#!/usr/bin/env python3
"""
Convenience CLI wrapper for `mapping_validation`.

Usage examples (from repo root or `bcp/`):

    python bcp/validate_mapping_cli.py \
        --mapping bcp/docs/NVUS2024101701-29-mapping.csv \
        --provider novogene \
        --data raw \
        --assay scale \
        --sif bcp/docs/SIF.csv

    python bcp/validate_mapping_cli.py \
        --mapping bcp/docs/NVUS2024101701-29-mapping.csv \
        --provider novogene \
        --data raw \
        --assay 10x

This simply forwards all arguments to `mapping_validation.main()`.
"""

from mapping_validation import main


if __name__ == "__main__":
    main()

