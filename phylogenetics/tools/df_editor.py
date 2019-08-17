__description__ = \
"""
Tools for editing sequences in phylopandas dataframes.
"""
__author__ = "Michael J. Harms (harmsm@gmail.com)"
__date__ = "2019-08-16"

import phylopandas as phy
import pandas as pd
import numpy as np

import random, string

def _check_uid(df,uid):
    """
    Make sure uid is in dataframe and is unique.
    """

    num_hits = np.sum(df.uid == uid)
    if num_hits == 0:
        err = "uid '{}' not in data frame.\n".format(uid)
        raise ValueError(err)
    elif num_hits > 1:
        err = "uid '{}' found multiple times in data frame.\n".format(uid)
        raise ValueError(err)
    else:
        pass


def trim_sequence(df,uid,start=0,end=None):
    """
    Trim a sequence at specified indexes.

    uid: unique phylopandas identifier for sequence.
    start: slice start
    end: slice end

    returns: copy of dataframe with trimmed sequence replacing current
             sequence
    """

    _check_uid(df,uid)

    new_df = df.copy()

    new_sequence = new_df[new_df.uid == uid].sequence.tolist()[0][start:end]
    i = new_df[new_df.uid == uid].index[0]
    new_df.at[i,"sequence"] = new_sequence

    return new_df

def split_sequence(df,uid,split_site):
    """
    Split a sequence in half, creating a new entry.  If a tree has already
    been loaded, this will create a polytomy at the new position.
    """

    _check_uid(df,uid)

    # Grab sequence to split
    sequence_to_split = df[df.uid == uid].sequence.tolist()[0]
    s1 = sequence_to_split[:split_site]
    s2 = sequence_to_split[split_site:]

    # Place to split frame
    i = np.where(df.uid == uid)[0][0]

    # Split frame
    df1 = df.iloc[:i].copy()
    df2 = df.iloc[np.array((i,i))].copy()
    df3 = df.iloc[(i+1):].copy()

    # Create new uid for new sequence
    uid_size = len(df.iloc[i].uid)
    new_uid = "".join([random.choice(string.ascii_letters)
                       for i in range(uid_size)])
    # Update entries with splits and new id
    df2 = df2.reset_index(drop=True)
    df2.at[0,"sequence"] = s1
    df2.at[1,"sequence"] = s2
    df2.at[1,"uid"] = new_uid

    # Concatenate and reindex
    new_df = pd.concat([df1,df2,df3],ignore_index=True,axis=0)
    new_df = new_df.reset_index(drop=True)

    return new_df
