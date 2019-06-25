import datetime as dt
from functools import partial
import numpy as np
import pandas as pd
from scipy import optimize
from scipy.special import expit
from scipy.special import logit

import utils


def distance_joint_distrib_three_periods(
    upsilon,
    upsilon_pseudo_inverse,
    joint_distrib_Y_three_periods,
    joint_distrib_Y_two_periods,
    diagonal_pr_Y_given_S,
):

    return np.linalg.norm(
        joint_distrib_Y_three_periods
        - np.dot(
            np.dot(np.dot(upsilon, diagonal_pr_Y_given_S), upsilon_pseudo_inverse),
            joint_distrib_Y_two_periods,
        ),
        ord="fro",
    )


def distance_joint_distrib_two_periods(upsilon, joint_distrib_Y, joint_distrib_truth):

    return np.linalg.norm(
        joint_distrib_Y
        - np.dot(np.dot(upsilon, joint_distrib_truth), np.transpose(upsilon)),
        ord="fro",
    )


def constraint_function(
    x,
    joint_distrib_Y_list,
    joint_distrib_Y_dict_three_periods,
    Y_sorted,
    S_sorted,
    third_period_land_uses,
):
    initial_distribution, pr_transition_list, upsilon = expand_objfn_x(
        x, Y_sorted, S_sorted, n_years=len(joint_distrib_Y_list) + 1
    )
    initial_distribution_sum = np.sum(initial_distribution)
    pr_transition_row_sums = np.hstack(
        [pr_transition.sum(axis=1).flatten() for pr_transition in pr_transition_list]
    )
    upsilon_column_sums = np.sum(upsilon, axis=0)
    return np.hstack(
        [
            initial_distribution_sum - 1.0,
            pr_transition_row_sums - 1.0,
            upsilon_column_sums - 1.0,
        ]
    )


def get_pr_transition_hat_list_as_dataframes(pr_transition_hat_list, S_sorted):

    for index in range(len(pr_transition_hat_list)):

        pr_transition_hat = pd.DataFrame(pr_transition_hat_list[index])
        pr_transition_hat.index = S_sorted
        pr_transition_hat.columns = S_sorted
        pr_transition_hat_list[
            index
        ] = pr_transition_hat  # Initially an array, now a dataframe

    return pr_transition_hat_list


def get_minimum_distance_parameter_estimates(
    dataset,
    df,
    Y_sorted,
    S_sorted,
    years,
    x_initial,
    third_period_land_uses,
    max_iterations,
):

    print("running minimum distance estimation, time is {}".format(dt.datetime.now()))
    print(
        " third period land uses used in estimation ({}): {}".format(
            len(third_period_land_uses), ", ".join(third_period_land_uses)
        )
    )

    joint_distrib_Y_list = list()
    for year in years[0:-1]:
        crosstab = utils.get_crosstab(dataset, df, Y_sorted, year_from=year)
        assert crosstab.values.sum() == len(
            df
        )  # Crosstab values sum up to number of raster cells
        joint_distrib_Y_list.append(crosstab / len(df))  # Sums to 1.0

    assert np.all(
        np.isclose(
            joint_distrib_Y_list[0].sum(axis=1), joint_distrib_Y_list[1].sum(axis=0)
        )
    )  # Sanity check: same marginal distribution

    ## Similar to above, but now conditional on third period land use
    joint_distrib_dict_third_period = {}
    for land_use in third_period_land_uses:
        joint_distrib_Y_list_third_period = list()
        for year in years[0:-2]:
            crosstab = utils.get_crosstab(
                dataset, df, Y_sorted, year_from=year, third_period_land_use=land_use
            )
            assert crosstab.values.sum() < len(
                df
            )  # Smaller because conditioning on third period land use
            joint_distrib_Y_list_third_period.append(
                crosstab / len(df)
            )  # Sums to < 1.0

        joint_distrib_dict_third_period[land_use] = joint_distrib_Y_list_third_period

    wrapper_args = (
        joint_distrib_Y_list,
        joint_distrib_dict_third_period,  # Conditional on third period land use (dict keys)
        Y_sorted,
        S_sorted,
        third_period_land_uses,
    )
    print(" starting minimizer, time is {}".format(dt.datetime.now()))

    objfn_value_initial = objective_fn_minimum_distance_wrapper(
        x_initial, *wrapper_args
    )

    ## TODO Positive directional derivative for linesearch, try writing analytical gradient?
    minimizer_result = optimize.minimize(
        fun=objective_fn_minimum_distance_wrapper,
        x0=logit(x_initial),
        method="SLSQP",  # Sequential Least SQuares Programming
        constraints={
            "type": "eq",  # TODO 'jac'
            "fun": constraint_function,  # TODO Simple gradient
            "args": wrapper_args,
        },
        options={
            "maxiter": max_iterations,
            "eps": 1e-06,  # TODO Use defaults?
            "ftol": 1e-06,
            "disp": True,
        },
        args=wrapper_args,
    )  # Careful with args order
    ## TODO Save minimizer results, especially exit code (make sure it's zero)
    print(" done running minimizer, time is {}".format(dt.datetime.now()))

    objfn_value_final = objective_fn_minimum_distance_wrapper(
        minimizer_result["x"], *wrapper_args
    )
    print(
        " objective function values: initial={}, final={}".format(
            objfn_value_initial, objfn_value_final
        )
    )

    estimated_parameters = expand_objfn_x(
        minimizer_result["x"], Y_sorted, S_sorted, n_years=len(years)
    )
    initial_distribution_hat, pr_transition_hat_list, upsilon_hat = estimated_parameters
    pr_transition_hat_list = get_pr_transition_hat_list_as_dataframes(
        pr_transition_hat_list, S_sorted
    )
    return initial_distribution_hat, pr_transition_hat_list, upsilon_hat


def expand_objfn_x(x, Y_sorted, S_sorted, n_years, logistic=True):

    assert n_years > 1
    assert x.size == (
        len(S_sorted)
        + (n_years - 1) * (len(S_sorted) ** 2)  # Initial distribution
        + len(Y_sorted) * len(S_sorted)  # Transitions
    )  # Upsilon

    if logistic:

        # Note: expit is the logistic function expit(x) = 1/(1 + exp(-x))
        # Map from (-Inf, +Inf) to [0, 1]
        x = expit(x)

    initial_distribution = np.array((x[0 : len(S_sorted)]))
    pr_transition_list = list()
    for i in range(n_years - 1):
        index_start = len(S_sorted) + i * (len(S_sorted) ** 2)
        index_stop = len(S_sorted) + (i + 1) * (len(S_sorted) ** 2)
        pr_transition_list.append(
            np.array(x[index_start:index_stop]).reshape((len(S_sorted), len(S_sorted)))
        )

    upsilon = pd.DataFrame(
        np.array(x[-(len(S_sorted) * len(Y_sorted)) :]).reshape(
            (len(Y_sorted), len(S_sorted))
        )
    )
    upsilon.index = Y_sorted
    upsilon.columns = S_sorted
    return initial_distribution, pr_transition_list, upsilon


def get_list_diagonal_pr_Y_given_S(
    upsilon, joint_distrib_S_list, Y_sorted, S_sorted, third_period_land_use
):
    upsilon_third_period_land_use_row = np.expand_dims(
        upsilon.loc[third_period_land_use].values, axis=1
    )
    assert upsilon_third_period_land_use_row.size == len(S_sorted)
    upsilon_third_period_land_use = np.tile(
        upsilon_third_period_land_use_row, len(S_sorted)
    )
    assert np.all(
        upsilon_third_period_land_use[0, 0] == upsilon_third_period_land_use[0]
    )  # True land use along rows
    list_diagonal_pr_Y_given_S = list()
    for joint_distrib_truth in joint_distrib_S_list[1:]:
        assert (
            upsilon_third_period_land_use.shape == joint_distrib_truth.shape
        )  # |S|-by-|S|
        # if np.any(np.isclose(np.sum(joint_distrib_truth, axis=0), 0.0)):
        #     pdb.set_trace()  # TODO Lower bounds slightly above 0?
        ## TODO Saw np.sum(joint_distrib_truth) == 1.0000000000000053 > 1.0, slightly outside constraints?
        pr_observation_given_truth_numerator = np.sum(
            np.multiply(upsilon_third_period_land_use, joint_distrib_truth), axis=0
        )  # Element-wise
        pr_observation_given_truth_denominator = np.sum(joint_distrib_truth, axis=0)
        pr_observation_given_truth = np.zeros_like(
            pr_observation_given_truth_denominator
        )  # Pr[ Y_{t+2} = y | S_{t+1} ]
        denominator_nonzero = pr_observation_given_truth_denominator > 0.0
        pr_observation_given_truth[denominator_nonzero] = (
            pr_observation_given_truth_numerator[denominator_nonzero]
            / pr_observation_given_truth_denominator[denominator_nonzero]
        )
        diagonal_pr_Y_given_S = np.diag(
            pr_observation_given_truth
        )  # Zeros on off-diagonal
        # if np.any(diagonal_pr_Y_given_S > 1.0) or np.any(diagonal_pr_Y_given_S < 0.0):
        #     pdb.set_trace()  # TODO

        assert np.all(diagonal_pr_Y_given_S >= 0.0)
        # assert np.all(diagonal_pr_Y_given_S <= 1.0)  # TODO Failed with entry == 1.00000001
        list_diagonal_pr_Y_given_S.append(
            diagonal_pr_Y_given_S
        )  # Careful with indexing

    return list_diagonal_pr_Y_given_S


def get_joint_distrib_S_list(initial_distribution, pr_transition_list):
    joint_distrib_S_list = list()
    marginal_distribution = initial_distribution
    for index in range(len(pr_transition_list)):
        marginal_distribution_tiled = np.tile(
            marginal_distribution, (marginal_distribution.size, 1)
        )
        joint_distrib_S = marginal_distribution_tiled * pr_transition_list[index].T
        joint_distrib_S_list.append(joint_distrib_S)
        marginal_distribution = np.sum(joint_distrib_S, axis=1)

    return joint_distrib_S_list


def objective_fn_minimum_distance_wrapper(
    x,
    joint_distrib_Y_list,
    joint_distrib_Y_dict_three_periods,
    Y_sorted,
    S_sorted,
    third_period_land_uses,
    weighting_matrix=None,
):
    assert x.shape[0] == x.size  # Flattened version of (joint_distrib_S_list, upsilon)
    initial_distribution, pr_transition_list, upsilon = expand_objfn_x(
        x, Y_sorted, S_sorted, n_years=len(joint_distrib_Y_list) + 1
    )
    joint_distrib_S_list = get_joint_distrib_S_list(
        initial_distribution, pr_transition_list
    )
    return objective_fn_minimum_distance(
        joint_distrib_S_list=joint_distrib_S_list,
        upsilon=upsilon,
        joint_distrib_Y_list=joint_distrib_Y_list,
        joint_distrib_Y_dict_three_periods=joint_distrib_Y_dict_three_periods,
        Y_sorted=Y_sorted,
        S_sorted=S_sorted,
        third_period_land_uses=third_period_land_uses,
    )


def objective_fn_minimum_distance(
    joint_distrib_S_list,
    upsilon,
    joint_distrib_Y_list,
    joint_distrib_Y_dict_three_periods,
    Y_sorted,
    S_sorted,
    third_period_land_uses,
    weighting_matrix=None,
):
    assert len(joint_distrib_Y_list) == len(joint_distrib_S_list)
    list_of_distances_two_periods = list()
    for i in range(len(joint_distrib_S_list)):
        distance = distance_joint_distrib_two_periods(
            upsilon, joint_distrib_Y_list[i], joint_distrib_S_list[i]
        )
        list_of_distances_two_periods.append(distance)

    list_of_distances_three_periods = list()
    upsilon_pseudo_inverse = np.linalg.pinv(upsilon)  # TODO Careful, this can fail
    for third_period_land_use in third_period_land_uses:
        joint_distrib_Y_list_three_periods = joint_distrib_Y_dict_three_periods[
            third_period_land_use
        ]
        list_of_diagonal_pr_Y_given_S = get_list_diagonal_pr_Y_given_S(
            upsilon, joint_distrib_S_list, Y_sorted, S_sorted, third_period_land_use
        )
        assert len(list_of_diagonal_pr_Y_given_S) == len(
            joint_distrib_Y_list_three_periods
        )
        for i in range(len(joint_distrib_Y_list_three_periods)):
            distance = distance_joint_distrib_three_periods(
                upsilon,
                upsilon_pseudo_inverse,
                joint_distrib_Y_list_three_periods[i],
                joint_distrib_Y_list[i],
                list_of_diagonal_pr_Y_given_S[i],
            )  # Careful with index
            list_of_distances_three_periods.append(distance)

    distance = np.sum(
        np.array(list_of_distances_two_periods + list_of_distances_three_periods) ** 2
    )  # TODO weighting matrix, this assumes identity
    return distance


def get_em_parameter_estimates(
    dataset,
    df,
    Y_sorted,
    S_sorted,
    years,
    upsilon_hat,
    initial_distribution_hat,
    pr_transition_hat_list,
    max_iterations=10,
):

    print("running EM estimation, time is {}".format(dt.datetime.now()))
    observation_cols = [utils.get_observation_colname(year) for year in years]

    for iteration in range(max_iterations):
        baum_welch_arrays = get_baum_welch_arrays(
            df,
            S_sorted,
            years,
            observation_cols,
            upsilon_hat,
            initial_distribution_hat,
            pr_transition_hat_list,
        )

        initial_distribution_hat_previous = initial_distribution_hat  # For distance
        initial_distribution_hat = np.mean(baum_welch_arrays["pi"][:, :, 0], axis=0)
        assert utils.is_valid_initial_distribution(initial_distribution_hat, S_sorted)

        pr_transition_hat_list_previous = pr_transition_hat_list
        pr_transition_hat_list = []
        for t in range(len(years) - 1):
            joint_distribution_mean = np.mean(
                baum_welch_arrays["pi_transition"][:, t], axis=0
            )
            marginal_distribution_mean = np.mean(
                baum_welch_arrays["pi"][:, :, t], axis=0
            )
            marginal_distribution_mean = np.tile(
                np.expand_dims(marginal_distribution_mean, axis=1),
                len(marginal_distribution_mean),
            )
            assert marginal_distribution_mean.shape == joint_distribution_mean.shape
            pr_transition_hat = joint_distribution_mean / marginal_distribution_mean
            pr_transition_hat_list.append(pr_transition_hat)

        assert utils.is_valid_pr_transition_list(
            pr_transition_hat_list, S_sorted, years
        )

        upsilon_hat_previous = upsilon_hat.copy()
        for s_index in range(len(S_sorted)):
            for y_index in range(len(Y_sorted)):
                numerator = np.sum(
                    baum_welch_arrays["pi"][:, s_index, :]
                    * (df[observation_cols].values == Y_sorted[y_index])
                )
                denominator = np.sum(baum_welch_arrays["pi"][:, s_index, :])
                upsilon_hat.iloc[y_index, s_index] = numerator / denominator

        assert utils.is_valid_upsilon(upsilon_hat, Y_sorted, S_sorted)

        distance_upsilon = np.max(
            np.abs(upsilon_hat.values.flatten() - upsilon_hat_previous.values.flatten())
        )
        distance_pr_transition = np.max(
            np.abs(
                np.array(pr_transition_hat_list).flatten()
                - np.array(pr_transition_hat_list_previous).flatten()
            )
        )
        distance_initial = np.max(
            np.abs(initial_distribution_hat - initial_distribution_hat_previous)
        )
        format_string = (
            " EM iteration {} distances: upsilon {}, transitions {}, initial distrib {}"
        )
        print(
            format_string.format(
                iteration,
                np.round(distance_upsilon, 5),
                np.round(distance_pr_transition, 5),
                np.round(distance_initial, 5),
            )
        )  # TODO Early stopping if distance < epsilon
        print(" time is {}".format(dt.datetime.now()))

    pr_transition_hat_list = get_pr_transition_hat_list_as_dataframes(
        pr_transition_hat_list, S_sorted
    )
    return initial_distribution_hat, pr_transition_hat_list, upsilon_hat


def get_baum_welch_arrays(
    df,
    S_sorted,
    years,
    observation_cols,
    upsilon,
    initial_distribution,
    pr_transition_list,
):
    map_observations_to_baum_welch = {}
    log_likelihood = np.zeros((len(df),))
    pi = np.zeros(
        (len(df), len(S_sorted), len(years))
    )  # Posteriors over hidden state at time t
    pi_transition = np.zeros(
        (len(df), len(years) - 1, len(S_sorted), len(S_sorted))
    )  # Posteriors over joint distribution of hidden state at (t, t+1)
    for i in range(len(df)):
        observations = tuple(
            df.iloc[i][observation_cols].tolist()
        )  # Tuples can be used as dictionary keys
        if not observations in list(map_observations_to_baum_welch.keys()):
            baum_welch = get_baum_welch(
                observations=observations,
                upsilon=upsilon,
                initial_distribution=initial_distribution,
                pr_transition_list=pr_transition_list,
            )
            map_observations_to_baum_welch.update({observations: baum_welch})

        else:
            baum_welch = map_observations_to_baum_welch[observations]

        log_likelihood[i] = baum_welch["log_likelihood"]
        pi[i] = baum_welch["pi"]
        for t in range(len(years) - 1):
            pi_transition[i, t] = baum_welch["pi_transition_list"][t]

    return {"log_likelihood": log_likelihood, "pi": pi, "pi_transition": pi_transition}


def get_baum_welch(observations, upsilon, initial_distribution, pr_transition_list):
    ## Written following Ramon van Handel's HMM notes, page 40, algorithm 3.2
    ## https://www.princeton.edu/~rvan/orf557/hmm080728.pdf
    ## TODO Same as viterbi: save a dict mapping observations to posteriors, in case obs seen multiple times
    likelihoods = np.zeros((len(observations)))
    likelihoods[0] = np.sum(initial_distribution * upsilon.loc[observations[0]])
    ## Probabilities over hidden state s_t conditional on {y_0, y_1, ... , y_t}, as opposed to full history
    pi_contemporaneous = np.zeros((len(initial_distribution), len(observations)))
    pi_contemporaneous[:, 0] = (
        upsilon.loc[observations[0]] * initial_distribution / likelihoods[0]
    )
    ## Forward loop
    for t in range(1, len(observations)):
        upsilon_row = upsilon.loc[observations[t]]
        P_transpose = pr_transition_list[t - 1].T
        pi_tilde = upsilon_row * np.dot(P_transpose, pi_contemporaneous[:, t - 1])
        likelihoods[t] = np.sum(pi_tilde)
        pi_contemporaneous[:, t] = pi_tilde / likelihoods[t]

    beta = np.zeros_like(
        pi_contemporaneous
    )  # Number of hidden states by observation length
    beta[:, -1] = 1 / likelihoods[-1]
    ## Probabilities over hidden state s_t conditional on {y_0, y_1, ... y_t, ... , y_T}, i.e. full history
    pi = np.zeros_like(pi_contemporaneous)
    pi[:, -1] = pi_contemporaneous[
        :, -1
    ]  # History up until last period is, by definition, the full history
    pi_transition_list = (
        []
    )  # List of posterior probabilities over hidden state transitions
    ## Backward loop
    for t in reversed(list(range(0, len(observations) - 1))):
        upsilon_diag = np.diag(
            upsilon.loc[observations[t + 1]]
        )  # TODO Index error, +1?
        pi_diag = np.diag(pi_contemporaneous[:, t])
        beta_diag = np.diag(beta[:, t + 1])
        beta[:, t] = (
            np.dot(np.dot(pr_transition_list[t], upsilon_diag), beta[:, t + 1])
            / likelihoods[t]
        )
        pi_transition = np.dot(
            np.dot(np.dot(pi_diag, pr_transition_list[t]), upsilon_diag), beta_diag
        )
        pi_transition_list.append(
            pi_transition
        )  # Backwards, need to reverse after loop
        pi[:, t] = np.sum(pi_transition, axis=1)

    pi_transition_list.reverse()  # In-place?
    log_likelihood = np.sum(np.log(likelihoods))
    assert np.all(np.isclose(pi.sum(axis=0), 1.0))
    return {
        "log_likelihood": log_likelihood,
        "pi": pi,
        "pi_transition_list": pi_transition_list,
    }


def get_viterbi(observations, upsilon, initial_distribution, pr_transition_list):

    n_time_periods = len(pr_transition_list) + 1
    assert len(observations) == n_time_periods

    # Value function, log likelihood
    value = np.full((n_time_periods, len(initial_distribution)), -np.inf)

    # Note: avoid RuntimeWarning: divide by zero encountered in log
    indices_nonzero = np.where(
        np.logical_and(initial_distribution > 0.0, upsilon.loc[observations[0]] > 0.0)
    )[0]
    format_string = "Value function -Inf in entire first row, observations = {}"
    assert indices_nonzero.size > 0, format_string.format(
        ", ".join(observations.tolist())
    )
    value[0, indices_nonzero] = np.log(initial_distribution[indices_nonzero]) + np.log(
        upsilon.loc[observations[0]][indices_nonzero]
    )  # Upsilon needs row & colnames
    argmax_array = np.zeros_like(
        value, dtype=int
    )  # Denoted b in van Handel's PDF, Algorithm 3.4 on page 46

    for time_index in range(n_time_periods)[1:]:
        for i in range(len(initial_distribution)):
            log_transition_probs = np.full((len(initial_distribution),), -np.inf)
            indices_nonzero = np.where(
                pr_transition_list[time_index - 1].iloc[:, i] > 0.0
            )[0]
            if indices_nonzero.size == 0:
                continue  # Don't take logs of transition probabilities if they're all zero

            transition_probs_nonzero = pr_transition_list[time_index - 1].iloc[:, i][
                indices_nonzero
            ]
            log_transition_probs[indices_nonzero] = np.log(transition_probs_nonzero)
            argmax_array[time_index, i] = np.nanargmax(
                np.array(value[time_index - 1] + log_transition_probs)
            )
            transition_prob = pr_transition_list[time_index - 1].iloc[
                argmax_array[time_index, i], i
            ]
            observation_prob = upsilon.loc[observations[time_index]][i]
            if transition_prob > 0.0 and observation_prob > 0.0:
                transition_probs_argmax = pr_transition_list[time_index - 1].iloc[
                    argmax_array[time_index, i], i
                ]
                value[time_index, i] = (
                    value[time_index - 1, argmax_array[time_index, i]]
                    + np.log(transition_probs_argmax)
                    + np.log(upsilon.loc[observations[time_index]][i])
                )

        format_string = (
            "Value function -Inf in entire row at index {}, observations = {}"
        )
        assert not np.all(np.isinf(value[time_index])), format_string.format(
            time_index, ", ".join(observations)
        )

    assert np.all(
        np.any(np.isfinite(value), axis=1)
    )  # All rows contain at least one finite value
    assert np.all(np.isfinite(argmax_array)) and not np.any(np.isnan(argmax_array))
    viterbi_path_indices = np.zeros((len(observations),), dtype=int)  # Will overwite
    viterbi_path_indices[-1] = np.argmax(value[-1])
    for time_index in range(n_time_periods)[1:]:
        viterbi_path_indices[-1 - time_index] = argmax_array[
            -time_index, viterbi_path_indices[-time_index]
        ]

    # Maximum likelihood sequence of ground truth landuses conditional on observations
    return viterbi_path_indices


def get_viterbi_df(
    df,
    S_sorted,
    years,
    observation_cols,
    upsilon,
    initial_distribution,
    pr_transition_list,
):

    # If we've already seen an observation sequence, can re-use its viterbi output
    map_observations_to_viterbi = {}

    viterbi_path_list = list()

    for i in range(len(df)):
        observations = tuple(
            df.iloc[i][observation_cols].tolist()
        )  # Tuples can be used as dictionary keys
        if observations in list(map_observations_to_viterbi.keys()):
            viterbi_path_list.append(map_observations_to_viterbi[observations])
        else:
            viterbi_land_uses = S_sorted[
                get_viterbi(
                    observations=observations,
                    upsilon=upsilon,
                    initial_distribution=initial_distribution,
                    pr_transition_list=pr_transition_list,
                )
            ]
            viterbi_path_list.append(viterbi_land_uses)
            map_observations_to_viterbi.update({observations: viterbi_land_uses})

    viterbi_paths_df = pd.DataFrame(np.array(viterbi_path_list))
    viterbi_paths_df.columns = ["viterbi_{}".format(x) for x in years]

    return viterbi_paths_df
