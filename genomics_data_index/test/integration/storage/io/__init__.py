from pathlib import Path

import pandas as pd
import vcfpy


def read_vcf_df(file: Path) -> pd.DataFrame:
    reader = vcfpy.Reader.from_path(path=str(file))
    df = pd.DataFrame([vars(r) for r in reader])
    df['TYPE'] = df['INFO'].apply(lambda x: x['TYPE'][0] if 'TYPE' in x else pd.NA)

    return df
