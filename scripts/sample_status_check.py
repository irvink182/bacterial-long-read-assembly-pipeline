import pandas as pd

df = pd.read_csv("pipeline_status.tsv", sep="\t")

def classify(sample_df):
    if (sample_df["status"] == "FAILED").any():
        return "FAILED"
    elif (sample_df["status"] != "OK").any():
        return "PARTIAL"
    else:
        return "OK"

final = (
    df.groupby("sample")
    .apply(classify)
    .reset_index(name="final_status")
)

final.to_csv("pipeline_final_status.tsv", sep="\t", index=False)