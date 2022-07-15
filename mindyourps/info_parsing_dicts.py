# File information
# ================

file_names_style = {
    # g1, genotype; g2, region
    'style1': r'(B3|D4|D8).*?(N450_MFNCR|N450|MFNCR)',
    # g1, dataset; g2, model settings (for BEAST runs)
    'style2': r'([^-]+)-([A-Z])\..*',
    # g1, genotype; g2, dataset; g3, region
    'style3': r'(?P<Genotype>B3|D4|D8)(?:_|-).*?'
              r'(?P<Dataset>GBR_GB|GBR_ROU_GB|ROU_GB|GB)(?:_|-)?.*?'
              r'(?P<Region>N450_MFNCR|N450|MFNCR)(?:_|-)?.*?'
              r'(?P<Settings>[A-Z])?',
    # tree summary files name format
    'style4': r'(?P<Genotype>B3|D4|D8)-(?P<Dataset>GBR_GB|GBR_ROU_GB|GB)-'
              r'(?P<Settings>[A-Z])',
    # stats files name format
    'style5': r'(?P<Dataset>GBR_GB|GBR_ROU_GB|GB)-(?P<Settings>[A-Z])-'
              r'(?P<Pivot_on>.*)-stats.csv'
}
