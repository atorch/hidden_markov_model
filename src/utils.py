import matplotlib

matplotlib.use(
    "Agg"
)  # See https://stackoverflow.com/questions/2801882/generating-a-png-with-matplotlib-when-display-is-undefined

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd


def get_crosstab(dataset, df, Y_sorted, year_from=2012, third_period_land_use=None):

    ## Careful with Y_sorted order: need crosstab to follow identical order in both rows and columns
    year_to = year_from + 1
    varname_from = "{}_{}".format(dataset, year_from)
    varname_to = "{}_{}".format(dataset, year_to)
    assert varname_from in df.columns and varname_to in df.columns
    if third_period_land_use is None:
        crosstab = pd.crosstab(
            df[varname_to], df[varname_from]
        )  # Careful: year_to along rows

    else:
        varname_third_period = "{}_{}".format(dataset, year_to + 1)
        df_subset = df.loc[df[varname_third_period] == third_period_land_use]
        if len(df_subset) == 0:
            ## Third period land use is never observed: crosstab is empty (all zeros)
            crosstab = pd.DataFrame(
                np.zeros((len(Y_sorted), len(Y_sorted))),
                index=Y_sorted,
                columns=Y_sorted,
            )
        else:
            crosstab = pd.crosstab(df_subset[varname_to], df_subset[varname_from])

    for land_use in Y_sorted:
        if not land_use in crosstab.columns:
            crosstab[land_use] = 0
        if not land_use in crosstab.index:
            crosstab.loc[land_use] = 0

    crosstab = (
        crosstab.sort_index()
    )  # Sort rows (so that they are in same order as Y_sorted)
    crosstab = crosstab.reindex_axis(
        sorted(crosstab.columns), axis=1
    )  # Sort columns (also same order as Y_sorted)

    assert crosstab.shape == (len(Y_sorted), len(Y_sorted))
    assert np.all(crosstab.index.tolist() == Y_sorted)
    assert np.all(crosstab.columns.tolist() == Y_sorted)

    return crosstab


def get_pr_transition_from_crosstab(crosstab):
    ## Careful, crosstab can contain empty (i.e. all zero) rows and/or columns
    ## Careful, crosstab has year_to along rows
    ## (transpose of my transition probability matrices, which have year_to along the column)
    crosstab_transpose = np.transpose(crosstab)
    row_sums = np.transpose(
        np.tile(np.sum(crosstab_transpose, axis=1), (crosstab_transpose.shape[1], 1))
    )
    pr_transition = crosstab_transpose / row_sums
    entire_row_is_nan = np.all(np.isnan(pr_transition), axis=1)
    assert np.all(
        np.logical_or(
            np.isclose(
                np.sum(pr_transition, axis=1), 1.0
            ),  # Sanity check: probabilities sum to 1
            entire_row_is_nan,
        )
    )  # ...or the entire row is NaN (land use not observed in given year)
    return (
        pr_transition
    )  # From along rows, to along columns -- careful, transpose of crosstab


def get_upsilon_diag_indices(Y_sorted, S_sorted):
    ## Indices for "diagonal" in non-square upsilon, i.e. (row, col) pairs where y == s
    diag_indices = (
        np.hstack([np.where(s == Y_sorted)[0] for s in S_sorted if s in Y_sorted]),
        np.arange(len(S_sorted)),
    )
    return diag_indices


def get_simple_upsilon(Y_sorted, S_sorted, diag_probabilities):

    ## Simple matrix of conditional observation probabilities (rows are classifications, columns are true land use)
    assert len(diag_probabilities) == len(S_sorted)
    assert len(S_sorted) <= len(Y_sorted)

    diag_probabilities = np.array(diag_probabilities)
    upsilon = np.zeros((len(Y_sorted), len(S_sorted)))
    diag_indices = get_upsilon_diag_indices(Y_sorted, S_sorted)
    upsilon[diag_indices] = diag_probabilities
    mask_off_diagonal = np.ones_like(upsilon, dtype=bool)
    mask_off_diagonal[diag_indices] = False
    off_diagonal_probabilities = (1.0 - diag_probabilities) / (len(Y_sorted) - 1)
    values_off_diagonal = (
        np.repeat(off_diagonal_probabilities, len(Y_sorted))
        .reshape(len(S_sorted), len(Y_sorted))
        .transpose()
    )
    upsilon[mask_off_diagonal] = values_off_diagonal[mask_off_diagonal]
    assert is_valid_upsilon(upsilon, Y_sorted, S_sorted)

    df_upsilon = pd.DataFrame(upsilon)
    df_upsilon.index = Y_sorted
    df_upsilon.columns = S_sorted

    return df_upsilon


def is_valid_years(years):
    return all(
        years[i] == years[i + 1] - 1 for i in range(len(years) - 1)
    )  # Sorted and increasing by 1


def is_valid_initial_distribution(initial_distribution, S_sorted):
    return (
        len(S_sorted) == len(initial_distribution)
        and np.all(initial_distribution >= 0.0)
        and np.isclose(  # Probabilities are not negative
            np.sum(initial_distribution), 1.0
        )
    )  # ...and they sum to 1.0


def is_valid_pr_transition_list(pr_transition_list, S_sorted, years):

    each_transition_matrix_is_valid = True

    for pr_transition in pr_transition_list:
        if (
            pr_transition.shape[0] != len(S_sorted)
            or pr_transition.shape[1] != len(S_sorted)
            or not np.all(np.isclose(np.sum(pr_transition, axis=1), 1.0))
        ):
            each_transition_matrix_is_valid = False

    return (
        len(pr_transition_list) == (len(years) - 1) and each_transition_matrix_is_valid
    )


def is_valid_upsilon(upsilon, Y_sorted, S_sorted):

    return (
        upsilon.shape == (len(Y_sorted), len(S_sorted))
        and np.all(upsilon >= 0.0)
        and np.all(np.isclose(np.sum(upsilon, axis=0), 1.0))
    )


def get_directory(
    dataset, years, aoi_width, aoi_height, aoi_x_min=None, aoi_y_min=None
):
    if "simulation" in dataset:
        ## No AOI coordinates for simulation, but width and height control sample size
        format_string = "{}_{}_to_{}_width_{}_height_{}"
        return format_string.format(
            dataset, min(years), max(years), aoi_width, aoi_height
        )
    else:
        format_string = "{}_{}_to_{}_x_min_{}_width_{}_y_min_{}_height_{}"
        return format_string.format(
            dataset, min(years), max(years), aoi_x_min, aoi_width, aoi_y_min, aoi_height
        )


def get_outfile(
    dataset,
    years,
    aoi_width,
    aoi_height,
    aoi_x_min,
    aoi_y_min,
    parent_directory,
    filename,
):
    subdirectory = get_directory(
        dataset, years, aoi_width, aoi_height, aoi_x_min, aoi_y_min
    )
    if not os.path.exists(parent_directory):
        os.makedirs(parent_directory)

    directory = os.path.join(parent_directory, subdirectory)
    if not os.path.exists(directory):
        os.makedirs(directory)

    outfile = os.path.join(directory, filename)
    return outfile


def plot_pr_transition(
    dataset,
    algorithm,
    years,
    pr_transition_hat_list,
    pr_transition_Y_list,
    pr_transition_truth_list,  # Only known for simulation
    land_use_from,
    land_use_to=None,
    aoi_x_min=None,
    aoi_width=None,
    aoi_y_min=None,
    aoi_height=None,
    linewidth=5.0,
    marker_size=14.0,
    text_size=18.0,
    plot_directory="./plots",
):
    assert len(years) == (len(pr_transition_hat_list) + 1)
    assert len(pr_transition_hat_list) == len(pr_transition_Y_list)

    color_map = (
        plt.cm.Dark2
    )  # See https://matplotlib.org/examples/color/colormaps_reference.html
    fig, ax = plt.subplots(figsize=(20, 10))
    ax.ticklabel_format(useOffset=False)  # Don't display years as e.g. '1 + 2.009e3'
    ax.set_yticks(np.arange(0.0, 1.01, 0.10), minor=False)
    ax.set_yticks(np.arange(0.05, 0.96, 0.10), minor=True)
    ax.grid(color="grey", linestyle="--", axis="y", which="major", alpha=0.75)
    ax.grid(color="grey", linestyle="--", axis="y", which="minor", alpha=0.5)

    if land_use_to is None:
        land_use_to = land_use_from  # E.g. forest-to-forest transition probabilities

    pr_transition_hat = [
        p.loc[land_use_from, land_use_to] for p in pr_transition_hat_list
    ]
    plt.plot(
        years[:-1],
        pr_transition_hat,
        color=color_map(0.0),
        linewidth=linewidth,
        marker="o",
        ms=marker_size,
    )
    ax.text(
        years[-1] - 0.80,
        pr_transition_hat[-1],
        "hmm estimates",
        color=color_map(0.0),
        verticalalignment="center",
        fontsize=text_size,
    )

    pr_transition_Y = [p.loc[land_use_from, land_use_to] for p in pr_transition_Y_list]
    plt.plot(
        years[:-1],
        pr_transition_Y,
        label="observations",
        color=color_map(0.2),
        linewidth=linewidth,
        marker="o",
        ms=marker_size,
    )
    ax.text(
        years[-1] - 0.80,
        pr_transition_Y[-1],
        "raster observations",
        color=color_map(0.2),
        verticalalignment="center",
        fontsize=text_size,
    )

    if pr_transition_truth_list is not None:
        assert (
            "simulation" in dataset
        )  # Sanity check, e.g. dataset == 'simulation_simple'
        pr_transition_truth = [
            p.loc[land_use_from, land_use_to] for p in pr_transition_truth_list
        ]
        plt.plot(
            years[:-1],
            pr_transition_truth,
            label="truth (simulation params)",
            color="grey",
            linewidth=linewidth,
            marker="o",
            ms=marker_size,
            linestyle="--",
        )
        offset = -0.02 + 0.04 * (pr_transition_truth[0] > pr_transition_hat[0])
        ax.text(
            years[0] + 0.20,
            pr_transition_truth[0] + offset,
            "truth (simulation params)",
            color="grey",
            verticalalignment="center",
            fontsize=text_size,
        )

    title_string = "{}: {}-to-{} transition probabilities".format(
        dataset, land_use_from, land_use_to
    )  # From t to t+1
    plt.title(title_string, fontsize=22)
    plt.ylabel("transition probability")
    plt.xlim(xmax=years[-1] + 0.50)  # Make room for ax.tex
    plt.xlim(xmin=years[0] - 0.25)
    plt.ylim(ymax=1.04)
    plt.ylim(ymin=min(min(pr_transition_hat), min(pr_transition_Y)) - 0.04)

    filename = "{}_pr_transition_{}_to_{}.png".format(
        algorithm.replace(" ", "_"), land_use_from, land_use_to
    )
    outfile = get_outfile(
        dataset,
        years,
        aoi_width,
        aoi_height,
        aoi_x_min,
        aoi_y_min,
        plot_directory,
        filename,
    )
    print(" saving {}".format(outfile))
    plt.savefig(outfile)
    plt.close("all")
    return
