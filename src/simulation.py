import datetime as dt
import numpy as np
import pandas as pd

from estimation import get_minimum_distance_parameter_estimates, get_viterbi_df
import utils

SIMULATION_CORN_SOY = "simulation_corn_soy"
SIMULATION_SIMPLE = "simulation_simple"
SIMULATION_TIME_VARYING = "simulation_time_varying"

SIMULATIONS = [SIMULATION_CORN_SOY, SIMULATION_SIMPLE, SIMULATION_TIME_VARYING]


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


def get_observation_df_and_Y_sorted(
    simulation_type, params, years, aoi_width, aoi_height
):

    initial_distribution = params["initial_distribution"]
    pr_transition_list = params["pr_transition_list"]
    upsilon = params["upsilon"]
    Y_sorted = params["Y_sorted"]
    S_sorted = params["S_sorted"]
    assert utils.is_valid_years(years)
    assert utils.is_valid_initial_distribution(initial_distribution, S_sorted)
    df = simulate_land_cover_df(
        simulation_type,
        Y_sorted=Y_sorted,
        S_sorted=S_sorted,
        initial_distribution=initial_distribution,
        pr_transition_list=pr_transition_list,
        upsilon=upsilon,
        years=years,
        aoi_width=aoi_width,
        aoi_height=aoi_height,
    )  # To increase sample size, increase width and height
    return df, Y_sorted


def print_classification_accuracy(simulation_type, years, df_with_viterbi):
    accuracy_dict_simulation, accuracy_dict_viterbi = {}, {}
    for year in years:
        col_truth = "true_land_use_{}".format(year)
        col_raster = "{}_{}".format(simulation_type, year)
        col_viterbi = "viterbi_{}".format(year)
        accuracy_dict_simulation[year] = np.mean(
            df_with_viterbi[col_truth] == df_with_viterbi[col_raster]
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
    print(" " + " ".join(accuracy_df.to_string().splitlines(True)))

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
    print(" " + np.array_str(initial_distribution.round(4)))
    print("true upsilon used in simulation:")
    print(" " + " ".join(df_upsilon.round(4).to_string().splitlines(True)))
    print("true transition probabilities used in simulation:")
    for pr_transition in pr_transition_list:
        print(" " + " ".join(pr_transition.round(4).to_string().splitlines(True)))

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
        format_string.format(
            max_abs_error_initial_distribution,
            max_abs_error_upsilon,
            max_abs_errors_pr_transition_string,
        )
    )
    return


def simulate_land_cover_single_point(
    simulation_type,
    Y_sorted,
    S_sorted,
    initial_distribution,
    pr_transition_list,
    upsilon,
    years=list(range(2010, 2017)),
    x=1,
    y=1,
):

    ## Simulation over time (at a single point)
    dict_for_series = {"aoi_x": x, "aoi_y": y}
    for year_idx, year in enumerate(years):
        land_cover_varname = "true_land_use_{}".format(year)
        observation_varname = "{}_{}".format(simulation_type, year)
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
        dict_for_series.update(
            {land_cover_varname: true_land_use, observation_varname: observation}
        )
        previous_land_cover = true_land_use

    row = pd.Series(dict_for_series)
    return row


def simulate_land_cover_df(
    simulation_type,
    Y_sorted,
    S_sorted,
    initial_distribution,
    pr_transition_list,
    upsilon,
    years=list(range(2010, 2017)),
    aoi_width=100,
    aoi_height=100,
):

    ## Simulation where true parameters are known, for sanity checking estimation code
    print("simulating raster data, time is {}".format(dt.datetime.now()))
    assert utils.is_valid_initial_distribution(initial_distribution, S_sorted)
    assert utils.is_valid_pr_transition_list(pr_transition_list, S_sorted, years)
    assert utils.is_valid_years(years)
    x_values = np.arange(0, aoi_width)
    y_values = np.arange(0, aoi_height)
    df = pd.DataFrame(
        np.transpose(
            [np.tile(x_values, len(y_values)), np.repeat(y_values, len(x_values))]
        )
    )  # Cartesian product
    df.columns = ["aoi_x", "aoi_y"]
    assert len(df) == len(x_values) * len(y_values)
    row_list = list()
    for index, row in df.iterrows():
        ## TODO This is somewhat slow, could be run in parallel (when model has no spatial dependence)
        simulation_row = simulate_land_cover_single_point(
            simulation_type,
            Y_sorted,
            S_sorted,
            initial_distribution,
            pr_transition_list,
            upsilon,
            years=years,
            x=row["aoi_x"],
            y=row["aoi_y"],
        )
        row_list.append(simulation_row)

    df = pd.DataFrame(row_list)

    return df


def main():
    print("hello")


if __name__ == "__main__":
    main()
