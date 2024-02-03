def infernal_scan_to_df(infernal):
    "DEPRECATED"

    # initiliaze an empty dataframe
    columns = [
        "hit_idx",
        "target_name",
        "t_accession",
        "query_name",
        "q_accession",
        "clan_name",
        "mdl",
        "mdl_from",
        "mdl_to",
        "seq_from",
        "seq_to",
        "strand",
        "trunc",
        "pass",
        "gc",
        "bias",
        "score",
        "E_value",
        "inc",
        "olp",
        "anyidx",
        "afrct1",
        "afrct2",
        "winidx",
        "wfrct1",
        "wfrct2",
        "description_of_target",
    ]

    df = pd.DataFrame(columns=columns)

    with open(infernal) as f:
        for line in f:
            if not line.startswith("#"):
                hit_list = line.strip().split(maxsplit=17)
                print(hit_list)
                d = {
                    "hit_idx": hit_list[0],
                    "target_name": hit_list[1],
                    "t_accession": hit_list[2],
                    "query_name": hit_list[3],
                    "q_accession": hit_list[4],
                    "clan_name": hit_list[5],
                    "mdl": hit_list[6],
                    "mdl_from": hit_list[7],
                    "mdl_to": hit_list[8],
                    "seq_from": hit_list[9],
                    "seq_to": hit_list[10],
                    "strand": hit_list[11],
                    "trunc": hit_list[12],
                    "pass": hit_list[13],
                    "gc": hit_list[14],
                    "bias": hit_list[15],
                    "score": hit_list[16],
                    "E_value": hit_list[17],
                    "inc": hit_list[18],
                    "olp": hit_list[19],
                    "anyidx": hit_list[20],
                    "afrct1": hit_list[21],
                    "afrct2": hit_list[22],
                    "winidx": hit_list[23],
                    "wfrct1": hit_list[24],
                    "wfrct2": hit_list[25],
                    "description_of_target": hit_list[26],
                }

                df = df.append(d, ignore_index=True)

                break

    df2 = df.drop(columns=["hit_idx", "anyidx", "afrct1", "winidx", "wfrct1", "wfrct2"])

    return df2


def infernal_search_to_df(infernal):
    # initiliaze an empty dataframe
    columns = [  # "hit_idx",
        "target_name",
        "t_accession",
        "query_name",
        "q_accession",
        # "clan_name",
        "mdl",
        "mdl_from",
        "mdl_to",
        "seq_from",
        "seq_to",
        "strand",
        "trunc",
        "pass",
        "gc",
        "bias",
        "score",
        "E_value",
        "inc",
        "description_of_target",
    ]

    df = pd.DataFrame(columns=columns)

    data = []
    with open(infernal) as f:
        for line in f:
            if not line.startswith("#"):
                hit_list = line.strip().split(maxsplit=17)

                d = {
                    "target_name": hit_list[0],
                    "t_accession": hit_list[1],
                    "query_name": hit_list[2],
                    "q_accession": hit_list[3],
                    "mdl": hit_list[4],
                    "mdl_from": int(hit_list[5]),
                    "mdl_to": int(hit_list[6]),
                    "seq_from": int(hit_list[7]),
                    "seq_to": int(hit_list[8]),
                    "strand": hit_list[9],
                    "trunc": hit_list[10],
                    "pass": hit_list[11],
                    "gc": float(hit_list[12]),
                    "bias": hit_list[13],
                    "score": float(hit_list[14]),
                    "E_value": float(hit_list[15]),
                    "inc": hit_list[16],
                    "description_of_target": hit_list[17],
                }

                # df = df.append(d, ignore_index = True)
                data.append(d)

    if len(data) > 0:
        df = pd.DataFrame.from_records(data)

    return df
