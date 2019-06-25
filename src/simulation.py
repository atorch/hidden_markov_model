import argparse
import datetime as dt
import numpy as np
import os
import pandas as pd
import pickle

from estimation import (
    get_em_parameter_estimates,
    get_minimum_distance_parameter_estimates,
    get_viterbi_df,
)
import utils

SIMULATION_CORN_SOY = "simulation_corn_soy"
SIMULATION_SIMPLE = "simulation_simple"
SIMULATION_TIME_VARYING = "simulation_time_varying"

SIMULATIONS = [SIMULATION_CORN_SOY, SIMULATION_SIMPLE, SIMULATION_TIME_VARYING]

OBSERVATIONS_TO_EXCLUDE_FROM_S = ["not_observed"]


def get_simulation_params(simulation_type, years):

    assert simulation_type in SIMULATIONS

    if simulation_type == SIMULATION_CORN_SOY:

        Y_sorted = np.sort(["corn", "forest", "not_observed", "pasture", "soy"])
        S_sorted = np.sort(["corn", "forest", "pasture", "soy"])
        initial_distribution = np.array([0.25, 0.25, 0.20, 0.30])
        upsilon = utils.get_simple_upsilon(
            Y_sorted, S_sorted, diag_probabilities=[0.85, 0.95, 0.80, 0.85]
        )
        pr_transition = pd.DataFrame(
            np.array(
                [
                    [0.05, 0.01, 0.04, 0.90],  # Transitions from corn
                    [0.01, 0.95, 0.03, 0.01],  # Transitions from forest
                    [0.02, 0.01, 0.95, 0.02],  # Transitions from pasture
                    [0.75, 0.01, 0.04, 0.20],
                ]
            )
        )  # Transitions from soy
        pr_transition.index = S_sorted
        pr_transition.columns = S_sorted
        pr_transition_list = [pr_transition] * (len(years) - 1)

    elif simulation_type == SIMULATION_SIMPLE:

        Y_sorted = np.sort(["crops", "forest", "pasture", "not_observed"])
        S_sorted = np.sort(["crops", "forest", "pasture"])
        initial_distribution = np.array([0.10, 0.60, 0.30])
        upsilon = utils.get_simple_upsilon(
            Y_sorted, S_sorted, diag_probabilities=[0.85, 0.90, 0.95]
        )
        upsilon.iloc[2, 0] = upsilon.iloc[2, 0] + upsilon.iloc[3:, 0].sum()
        upsilon.iloc[3:, 0] = 0.0
        assert utils.is_valid_upsilon(upsilon, Y_sorted, S_sorted)
        pr_transition = pd.DataFrame(
            np.array([[0.90, 0.01, 0.09], [0.01, 0.95, 0.04], [0.10, 0.05, 0.85]])
        )
        pr_transition.index = S_sorted
        pr_transition.columns = S_sorted
        pr_transition_list = [pr_transition] * (len(years) - 1)

    if simulation_type == SIMULATION_TIME_VARYING:

        Y_sorted = np.sort(["forest", "not_observed", "pasture", "soy"])
        S_sorted = np.sort(["forest", "pasture", "soy"])
        initial_distribution = np.array([0.50, 0.20, 0.30])
        upsilon = utils.get_simple_upsilon(
            Y_sorted, S_sorted, diag_probabilities=[0.95, 0.80, 0.70]
        )
        ## Time varying transition probabilities: one matrix in even years, another matrix in odd years
        pr_transition_even = pd.DataFrame(
            np.array(
                [
                    [0.96, 0.03, 0.01],  # Transitions from forest
                    [0.01, 0.95, 0.04],  # Transitions from pasture
                    [0.01, 0.04, 0.95],
                ]
            )
        )  # Transitions from soy
        pr_transition_even.index = S_sorted
        pr_transition_even.columns = S_sorted
        pr_transition_odd = pd.DataFrame(
            np.array(
                [
                    [0.90, 0.09, 0.01],  # Transitions from forest
                    [0.01, 0.98, 0.01],  # Transitions from pasture
                    [0.01, 0.09, 0.90],
                ]
            )
        )  # Transitions from soy
        pr_transition_odd.index = S_sorted
        pr_transition_odd.columns = S_sorted
        pr_transition_list = [
            pr_transition_even if year % 2 == 0 else pr_transition_odd
            for year in years[0:-1]
        ]

    simulation_params = {
        "Y_sorted": Y_sorted,
        "S_sorted": S_sorted,
        "initial_distribution": initial_distribution,
        "pr_transition_list": pr_transition_list,
        "upsilon": upsilon,
    }

    return simulation_params


def get_observation_df_Y_sorted(simulation_type, params, years):

    initial_distribution = params["initial_distribution"]
    pr_transition_list = params["pr_transition_list"]
    upsilon = params["upsilon"]
    Y_sorted = params["Y_sorted"]
    S_sorted = params["S_sorted"]
    assert utils.is_valid_years(years)
    assert utils.is_valid_initial_distribution(initial_distribution, S_sorted)

    # Note: to increase sample size, increase width and height
    df = simulate_land_cover_df(
        simulation_type,
        Y_sorted=Y_sorted,
        S_sorted=S_sorted,
        initial_distribution=initial_distribution,
        pr_transition_list=pr_transition_list,
        upsilon=upsilon,
        years=years,
    )

    return df, Y_sorted


def print_classification_accuracy(simulation_type, years, df_with_viterbi):

    accuracy_dict_simulation, accuracy_dict_viterbi = {}, {}

    for year in years:

        col_truth = "true_land_use_{}".format(year)
        col_observation = utils.get_observation_colname(year)
        col_viterbi = "viterbi_{}".format(year)
        accuracy_dict_simulation[year] = np.mean(
            df_with_viterbi[col_truth] == df_with_viterbi[col_observation]
        )
        accuracy_dict_viterbi[year] = np.mean(
            df_with_viterbi[col_truth] == df_with_viterbi[col_viterbi]
        )

    ## TODO This is accuracy by year, also print accuracy conditional on ground truth
    accuracy_df = pd.DataFrame(
        {
            "raster_accuracy": accuracy_dict_simulation,
            "viterbi_accuracy": accuracy_dict_viterbi,
        }
    )
    print("accuracy by year, for both original simulated raster and for viterbi:")
    print((" " + " ".join(accuracy_df.to_string().splitlines(True))))

    return


def print_parameter_estimate_accuracy(
    initial_distribution_hat,
    initial_distribution,
    pr_transition_list_hat,
    pr_transition_list,
    df_upsilon_hat,
    df_upsilon,
):

    print("true initial distribution used in simulation:")
    print((" " + np.array_str(initial_distribution.round(4))))
    print("true upsilon used in simulation:")
    print((" " + " ".join(df_upsilon.round(4).to_string().splitlines(True))))
    print("true transition probabilities used in simulation:")
    for pr_transition in pr_transition_list:
        print((" " + " ".join(pr_transition.round(4).to_string().splitlines(True))))

    max_abs_error_initial_distribution = np.max(
        np.abs(initial_distribution_hat - initial_distribution)
    )
    max_abs_error_upsilon = np.max(np.abs(df_upsilon_hat - df_upsilon).values)
    max_abs_errors_pr_transition = [
        np.max(np.abs(pr_transition_hat - pr_transition).values)
        for pr_transition_hat, pr_transition in zip(
            pr_transition_list_hat, pr_transition_list
        )
    ]
    max_abs_errors_pr_transition_string = ", ".join(
        ["{:0.4f}".format(x) for x in max_abs_errors_pr_transition]
    )
    format_string = "\n".join(
        [
            "maximum absolute errors in parameter estimates:",
            " initial distribution: {:0.4f}",
            " upsilon: {:0.4f}",
            " transition probs: {}",
        ]
    )
    print(
        (
            format_string.format(
                max_abs_error_initial_distribution,
                max_abs_error_upsilon,
                max_abs_errors_pr_transition_string,
            )
        )
    )
    return


def simulate_land_cover_single_pixel(
    simulation_type,
    Y_sorted,
    S_sorted,
    initial_distribution,
    pr_transition_list,
    upsilon,
    years,
    pixel_id,
):

    ## Simulation over time (at a single point/pixel)

    pixel = {"pixel_id": pixel_id}

    for year_idx, year in enumerate(years):

        land_cover_varname = "true_land_use_{}".format(year)
        observation_varname = utils.get_observation_colname(year)

        if year_idx == 0:
            true_land_use = np.random.choice(S_sorted, p=initial_distribution)
        else:
            probabilities = pr_transition_list[year_idx - 1].loc[
                previous_land_cover
            ]  # Markovian transitions
            true_land_use = np.random.choice(S_sorted, p=probabilities)

        observation = np.random.choice(
            Y_sorted, p=upsilon[true_land_use]
        )  # Columns of upsilon are Pr[Y | S]

        pixel.update(
            {land_cover_varname: true_land_use, observation_varname: observation}
        )
        previous_land_cover = true_land_use

    return pd.Series(pixel)


def simulate_land_cover_df(
    simulation_type,
    Y_sorted,
    S_sorted,
    initial_distribution,
    pr_transition_list,
    upsilon,
    years=list(range(2010, 2017)),
    n_pixels=1000,
):

    # Simulation where true parameters are known, for sanity checking estimation code

    print(("simulating raster data, time is {}".format(dt.datetime.now())))
    assert utils.is_valid_initial_distribution(initial_distribution, S_sorted)
    assert utils.is_valid_pr_transition_list(pr_transition_list, S_sorted, years)
    assert utils.is_valid_years(years)

    rows = list()

    for pixel_id in np.arange(n_pixels):

        # TODO This is slow, could be run in parallel (when model has no spatial dependence)
        row = simulate_land_cover_single_pixel(
            simulation_type,
            Y_sorted,
            S_sorted,
            initial_distribution,
            pr_transition_list,
            upsilon,
            years=years,
            pixel_id=pixel_id,
        )

        rows.append(row)

    return pd.DataFrame(rows)


def get_observation_df_Y_S(simulation_type, years, simulation_params):

    df, Y_sorted = get_observation_df_Y_sorted(
        simulation_type, simulation_params, years
    )

    print(
        f'done simulating data for {df["pixel_id"].nunique()} unique points across {len(years)} years'
    )

    print("set of possible observations Y (|Y| = {}):".format(len(Y_sorted)))
    print(" {}".format("\n ".join(Y_sorted)))

    marginal_list = []

    for year in years:

        observation_colname = utils.get_observation_colname(year)
        marginal = df[observation_colname].value_counts() / len(df)
        marginal_list.append(marginal)

    marginal_df = pd.DataFrame(
        marginal_list, index=years
    )  # Marginal distribution of Y by year
    marginal_df = marginal_df.fillna(value=0.0)

    # Note: the set of hidden states can be smaller than the set of observations
    S_sorted = [x for x in Y_sorted if x not in OBSERVATIONS_TO_EXCLUDE_FROM_S]

    return df, Y_sorted, S_sorted


def run_hmm_estimation(
    simulation_type, years, initial_pr_transition_diagonal, verbose=1
):

    assert len(years) >= 3, "Need at least 3 years of data to estimate HMM"

    simulation_params = get_simulation_params(simulation_type, years)

    df, Y_sorted, S_sorted = get_observation_df_Y_S(
        simulation_type, years, simulation_params=simulation_params
    )

    assert len(Y_sorted) >= len(S_sorted) > 1

    list_of_crosstabs = list()

    for year in years[0:-1]:
        list_of_crosstabs.append(
            utils.get_crosstab(simulation_type, df, Y_sorted, year_from=year)
        )

    crosstab_sum_equals_n_raster_cells = np.array(
        [crosstab.values.sum() for crosstab in list_of_crosstabs]
    ) == len(df)
    assert np.all(crosstab_sum_equals_n_raster_cells)
    pr_transition_Y_list = [
        utils.get_pr_transition_from_crosstab(crosstab)
        for crosstab in list_of_crosstabs
    ]

    ## Expect large off-diagonal raster transition probabilities, because of classification errors
    if verbose > 0:
        for index in range(len(years) - 1):
            print(
                "transition probabilities in raster observations, {} to {}:".format(
                    years[index], years[index + 1]
                )
            )
            print(
                " "
                + " ".join(
                    pr_transition_Y_list[index].round(4).to_string().splitlines(True)
                )
            )

    if not os.path.exists("./pickles"):
        os.makedirs("./pickles")

    pickle_outfile = utils.get_outfile(
        simulation_type,
        years,
        parent_directory="./pickles",
        filename="estimated_parameters.pickle",
    )

    if os.path.isfile(pickle_outfile):
        print(
            "{} already exists: will load pickle instead of re-running estimation".format(
                pickle_outfile
            )
        )
        estimates = pickle.load(open(pickle_outfile, "rb"))

    else:
        ## TODO Try eigenvalue decomp for initial values?
        initial_upsilon = utils.get_simple_upsilon(
            Y_sorted, S_sorted, diag_probabilities=0.95 * np.ones((len(S_sorted),))
        )

        # if 'crops_or_pasture' in Y_sorted and ('crops' in S_sorted or 'pasture' in S_sorted):
        #     pdb.set_trace()  # initial_upsilon.loc[np.logical_not(initial_upsilon.index == 'crops_or_pasture'), 'pasture']
        #     initial_upsilon  # TODO Set large initial Pr[ Y=crops_or_pasture | S=pasture ] ?

        ## TODO utils.get_simple_pr_transition
        off_diagonal = np.logical_not(np.eye(len(S_sorted), dtype=bool))
        initial_pr_transition = np.diag(
            initial_pr_transition_diagonal * np.ones((len(S_sorted),))
        )
        initial_pr_transition[off_diagonal] = (1 - initial_pr_transition_diagonal) / (
            len(S_sorted) - 1
        )
        initial_pr_transition_list = [initial_pr_transition] * (len(years) - 1)
        assert utils.is_valid_pr_transition_list(
            initial_pr_transition_list, S_sorted, years
        )
        initial_initial_distribution = np.ones((len(S_sorted),)) / float(
            len(S_sorted)
        )  # TODO initial_distrib_guess?

        ## TODO Also run EM estimates with constraints (e.g. Pr[Y = crops_or_pasture | S = old_forest] = 0), or transitions
        estimates_em = get_em_parameter_estimates(
            simulation_type,
            df,
            Y_sorted,
            S_sorted,
            years,
            initial_upsilon,
            initial_initial_distribution,
            initial_pr_transition_list,
        )

        x_initial = np.hstack(
            [initial_initial_distribution]
            + [pr_transition.flatten() for pr_transition in initial_pr_transition_list]
            + [initial_upsilon.values.flatten()]
        )
        third_period_land_uses = (
            Y_sorted
        )  # Optimization seems to work better when using more third period land uses
        estimates_min_dist = get_minimum_distance_parameter_estimates(
            simulation_type,
            df,
            Y_sorted,
            S_sorted,
            years,
            x_initial=x_initial,
            third_period_land_uses=third_period_land_uses,
            max_iterations=300,
        )  # TODO Argument
        print(
            "saving estimated parameters (both EM and minimum distance) to {}".format(
                pickle_outfile
            )
        )
        estimates = {"minimum distance": estimates_min_dist, "EM": estimates_em}
        pickle.dump(estimates, open(pickle_outfile, "wb"))

    for algorithm in list(estimates.keys()):
        initial_distribution_hat, pr_transition_hat_list, upsilon_hat = estimates[
            algorithm
        ]
        for index in range(len(years) - 1):
            print(
                "{} transition probabilities, {} to {}:".format(
                    algorithm, years[index], years[index + 1]
                )
            )
            print(
                " "
                + " ".join(
                    pr_transition_hat_list[index].round(4).to_string().splitlines(True)
                )
            )

        print("saving plots of {} transition probabilities".format(algorithm))
        pr_transition_truth_list = None

        pr_transition_truth_list = simulation_params[
            "pr_transition_list"
        ]  # For plot_pr_transition

        for land_use in S_sorted:
            utils.plot_pr_transition(
                simulation_type,
                algorithm,
                years,
                pr_transition_hat_list,
                pr_transition_Y_list,
                pr_transition_truth_list,
                land_use,
            )

        if simulation_type == SIMULATION_CORN_SOY:

            utils.plot_pr_transition(
                simulation_type,
                algorithm,
                years,
                pr_transition_hat_list,
                pr_transition_Y_list,
                pr_transition_truth_list,
                land_use_from="soy",
                land_use_to="corn",
            )

            utils.plot_pr_transition(
                simulation_type,
                algorithm,
                years,
                pr_transition_hat_list,
                pr_transition_Y_list,
                pr_transition_truth_list,
                land_use_from="corn",
                land_use_to="soy",
            )

        diag_indices = utils.get_upsilon_diag_indices(Y_sorted, S_sorted)
        diag_upsilon_string = ", ".join(
            np.round(upsilon_hat.values[diag_indices], 4).astype(str).tolist()
        )
        format_string = (
            "diagonals of observation probability matrix (upsilon) estimated by {}: {}"
        )
        print(format_string.format(algorithm, diag_upsilon_string))
        print("observation probability matrix estimated by {}:".format(algorithm))
        print(" " + " ".join(upsilon_hat.round(4).to_string().splitlines(True)))

        print_parameter_estimate_accuracy(
            initial_distribution_hat,
            simulation_params["initial_distribution"],
            pr_transition_hat_list,
            simulation_params["pr_transition_list"],
            upsilon_hat,
            simulation_params["upsilon"],
        )

        if not os.path.exists("./csv"):
            os.makedirs("./csv")

        filename = "df_with_viterbi_using_{}_estimates.csv".format(
            algorithm.replace(" ", "_")
        )
        viterbi_outfile = utils.get_outfile(
            simulation_type, years, parent_directory="./csv", filename=filename
        )

        if os.path.isfile(viterbi_outfile):
            print(
                "{} already exists: will load csv instead of re-running viterbi".format(
                    viterbi_outfile
                )
            )
            df_with_viterbi = pd.read_csv(viterbi_outfile)

        else:

            print(
                "running viterbi on raster data using {} estimates, time is {}".format(
                    algorithm, datetime.datetime.now()
                )
            )

            observation_cols = [utils.get_observation_colname(year) for year in years]

            viterbi_paths_df = get_viterbi_df(
                df,
                S_sorted,
                years,
                observation_cols,
                upsilon_hat,
                initial_distribution_hat,
                pr_transition_hat_list,
            )  # TODO Pass model parameters as a dictionary

            print(" done running viterbi, time is {}".format(datetime.datetime.now()))
            assert len(df) == len(viterbi_paths_df)
            df_with_viterbi = pd.concat([df, viterbi_paths_df], axis=1)
            assert len(df_with_viterbi) == len(df)
            viterbi_cols = ["viterbi_{}".format(x) for x in years]
            df_with_viterbi["raster_and_viterbi_disagree"] = np.any(
                df_with_viterbi[viterbi_cols].values
                != df_with_viterbi[observation_cols].values,
                axis=1,
            )
            print("saving {}".format(viterbi_outfile))
            df_with_viterbi.to_csv(viterbi_outfile, index=False, header=True)

            for year in years:

                observation_colname = utils.get_observation_colname(year)

                pr_viterbi_equals_raster = np.mean(
                    df_with_viterbi["viterbi_{}".format(year)]
                    == df_with_viterbi[observation_colname]
                )

                format_string = "viterbi agrees with {} raster in {} of rows"
                print(
                    format_string.format(year, pr_viterbi_equals_raster)
                )  # Often in [0.90, 0.95]

            print_classification_accuracy(
                simulation_type, years, df_with_viterbi
            )  # Uses viterbi


def main(simulation_type, year_start, year_end, initial_pr_transition_diagonal):

    # Note: code assumes years are sorted
    years = np.arange(year_start, year_end + 1)

    run_hmm_estimation(
        simulation_type=simulation_type,
        years=years,
        initial_pr_transition_diagonal=initial_pr_transition_diagonal,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--simulation_type", type=str, default=SIMULATION_SIMPLE, choices=SIMULATIONS
    )

    parser.add_argument(
        "--initial_pr_transition_diagonal",
        type=float,
        default=0.70,
        help="Initial value for transiton matrix diagonals (initial values for objective function)",
    )

    parser.add_argument(
        "--year_start", type=int, default=2011, help="First raster year, inclusive."
    )
    parser.add_argument(
        "--year_end", type=int, default=2016, help="Last raster year, inclusive."
    )

    args = parser.parse_args()
    assert args.year_start < args.year_end

    main(
        simulation_type=args.simulation_type,
        year_start=args.year_start,
        year_end=args.year_end,
        initial_pr_transition_diagonal=args.initial_pr_transition_diagonal,
    )
