"""Functions to parse and analyse MS dataset.

Functions:
    analyse_rawfile: Read excel, iso-correct and save as csv.
    parse_ms_data: Read MS excel and return pandas DataFrame.
    column_percentage: Calculate percentages of DataFrame.
"""
import pandas as pd
import picor as ic


def analyse_rawfile(
    infile,
    outfile=None,
    splitting_file=None,
    splitting_mapping=None,
    isotopologue_correction=True,
    resolution_correction=False,
):
    """Parse data, iso correct it and export as csv.

    :param infile: str or Path
        Path to input excel file
    :param outfile: str or path (default: None)
        Path to output file (should include ".csv" suffix)
    :param splitting_file: str or Path
        Path to excel file with label splitting data
    :param splitting_mapping: dict
        Sheet name of splitting_file with column label of MS data as values
    :return: pandas DataFrame
        (Iso-corrected) dataset
    """
    df_raw = parse_ms_data(infile)
    if isotopologue_correction:
        df = ic.calc_isotopologue_correction(
            df_raw,
            "K(ac)QLATK(ac)AAR",
            exclude_col=["Rep", "Exp"],
            resolution_correction=True,
        )
    else:
        df = df_raw
    if splitting_file and splitting_mapping:
        dfs_site = pd.read_excel(splitting_file, sheet_name=None)
        df = calc_site_fraction(df, dfs_site, splitting_mapping)
    elif splitting_file or splitting_mapping:
        raise ValueError("Both splitting_file and splitting_mapping have tobe set")
    if outfile:
        df.to_csv(outfile)
    return df


def parse_ms_data(infile):
    """Parse MS data and return pandas DataFrame.

    Parse excel file with MS data.
    Different labelled peptides in columns,
    Experiments and timepoints in rows
    :param infile: str or Path
        Filepath to excel file
    :returns: pandas DataFrame
    """
    columns = {
        "K(ac)QLATK(ac)AAR": "No label",
        "K(ac)QLATK(ac13C)AAR": "2C13",
        "K(ac)QLATK(ac*)AAR": "2C13 3H02",
        "K(ac13C)QLATK(ac13C)AAR": "4C13",
        "K(ac13C)QLATK(ac*)AAR": "4C13 3H02",
        "K(ac*)QLATK(ac*)AAR": "4C13 6H02",
    }
    df = pd.read_excel(infile, index_col=0)
    df.rename(columns=columns, inplace=True)
    return df


def calc_site_fraction(df_ms, dfs_site, mapping_sites, outfile=None):
    """Insert new site specific columns.

    Multiply site percentages with label column and insert new columns.
    Removes unspecific columns used for calculation.
    :param data: pandas DataFrame
        MS raw data with data points in rows and labels in columns
    :param dfs_site: pandas DataFrame
        Splitting data with columns with percentage starting with '%'
    :param splitting_mapping: dict
        Sheet name of splitting_file with column label of MS data as values
    :param outfile: str or path (default: None)
        Path to output file (should include ".csv" suffix)
    :return: pandas DataFrame
    """
    data = df_ms.copy()
    data = data.reset_index()
    for site, ms_col in mapping_sites.items():
        df_site = dfs_site[site]
        for col in df_site:
            if not col.startswith("% "):
                continue
            if pd.merge(data, df_site).isnull().sum().sum():
                raise ValueError("NaNs in either MS or site data.")
            site_col = " ".join([ms_col, col[2:]])
            data[site_col] = data[ms_col] * df_site[col] / 100
        data.drop(columns=ms_col, inplace=True)
    data = data.set_index("Time in h").sort_index(axis="columns")
    data = move_col_first(data, "No label")
    if outfile:
        print(f"Saving {outfile}")
        data.to_csv(outfile)
    return data


def column_percentage(df, usecols):
    """Calculate percentages across rows.

    Calculate percentages of columns so that each row adds up to 1.
    :param df: pandas DataFrame
        Data with column "Rep" and "Exp" column in addtion to usecols
    :param usecols: list of column names
        Columns to be used for calculation
    :returns: pandas DataFrame
    """
    df = df.loc[:, ["Exp", "Rep"] + usecols]
    df[usecols] = df[usecols].div(df[usecols].sum(axis=1), axis=0)  # .multiply(100)
    return df


def move_col_first(df, col):
    """Move col to first column"""
    df_order = df.copy()
    df_order.insert(0, col, df_order.pop(col))
    return df_order
