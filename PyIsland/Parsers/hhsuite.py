import pandas as pd


def parse_hhsearch_scores(scores, min_prob=20):
    """
    Converts the "--scores" output of HHSearch v1.5 to a dataframe
    stops parsing when probability < min_prob
    """
    i = 0
    data = []
    with open(scores) as infile:
        for line in infile:
            i += 1

            if i == 1:
                qseq = line.strip().split()[1]

            if i == 4:
                qlen = int(line.strip().split()[1])

            if i > 5:
                # print(line.strip().split())
                (
                    subject,
                    rel,
                    slen,
                    col,
                    logpva,
                    saass,
                    prob,
                    score,
                    logeval,
                ) = line.strip().split()

                if float(prob) < min_prob:
                    break

                data.append(
                    {
                        "qseq": qseq,
                        "qlen": qlen,
                        "subject": subject,
                        "rel": rel,
                        "slen": slen,
                        "col": int(col),
                        "logpva": logpva,
                        "saass": saass,
                        "prob": prob,
                        "score": score,
                        "logeval": logeval,
                    }
                )

    df = pd.DataFrame.from_records(data)

    if not df.empty:
        df2 = df.astype(
            {
                "qseq": "str",
                "qlen": "int64",
                "subject": "str",
                "slen": "int64",
                "col": "int64",
                "prob": "float",
            }
        )

        df2["qcov"] = df2["col"] / df2["qlen"]
        df2["scov"] = df2["col"] / df2["slen"]

    else:
        df2 = pd.DataFrame()

    return df2


def normalize_scores_by_selfhit(df):
    """
    Accepts a dataframe with columns qseq and subject
    Adds a column "score_norm" with the score divided by the self score
    returns a dataframe
    """

    num_selfhits = len(df.query("qseq == subject"))
    assert num_selfhits > 0, f"No self hits!"

    selfscore = float(df.query("qseq == subject").score.to_list()[0])

    df["score_norm"] = -np.log(df["score"].astype("float") / selfscore)

    return df



