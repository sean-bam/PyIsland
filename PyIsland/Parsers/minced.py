def minced_to_output(mincedgff, outcsv, outfasta):
    with open(outcsv, "w") as o:
        print(
            "contig",
            "start",
            "stop",
            "num_repeats",
            "crispr_id",
            "repeat",
            sep=",",
            file=o,
        )
        with open(mincedgff) as f:
            for line in f:
                if not line.startswith("#"):
                    # add the coordinates of the CRISPR array
                    (
                        contig,
                        drop1,
                        drop2,
                        start,
                        stop,
                        num_repeats,
                        drop4,
                        drop5,
                        crispr_id,
                    ) = line.strip().split()

                    crisprid, drop1, drop2, rpt = crispr_id.split(";")
                    crisprid = crisprid.split("=")[1]
                    rptseq = rpt.split("=")[1]

                    print(
                        contig,
                        start,
                        stop,
                        num_repeats,
                        crisprid,
                        rptseq,
                        sep=",",
                        file=o,
                    )

                    # make a fasta file of the repeats
                    with open(outfasta, "a") as o2:
                        crisprID_raw, drop2, drop3, repeat_seq_raw = crispr_id.split(
                            ";"
                        )
                        crisprID = crisprID_raw.split("=")[1]
                        repeat_seq = repeat_seq_raw.split("=")[1]
                        print(f">{contig}_{crisprID}", file=o2)
                        print(repeat_seq, file=o2)
