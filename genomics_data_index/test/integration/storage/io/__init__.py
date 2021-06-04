from pathlib import Path

import pandas as pd
import vcf


def read_vcf_df(file: Path) -> pd.DataFrame:
    reader = vcf.Reader(filename=str(file))
    df = pd.DataFrame([vars(r) for r in reader])
    df['TYPE'] = df['INFO'].apply(lambda x: x['TYPE'][0] if 'TYPE' in x else pd.NA)

    return df
