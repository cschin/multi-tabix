multi-tabix
================

`multi_tbx`: a simple tool for indexing VCF files and extract variant records for variant data stored in multiple VCF files.

If there are variant data scattered around multiple VCF files, the `multi_tbx` tool provides a way to make a meta-index by
scanning all tabix index files. With the meta-index, it can simplify the workflow to extract VCF records stored in multiple files.

For example, if we have multiple VCF files (assuming covering non-overlapped regions in a genome) and their associated `tbi` files:

```
❯ cat my_tbi_files
/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format/ukb23156_c10_b0_v1.vcf.gz.tbi
/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format/ukb23156_c10_b10_v1.vcf.gz.tbi
/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format/ukb23156_c10_b11_v1.vcf.gz.tbi
```

(We assume the vcf files are in the directory of the corresponded `tbi` files.)

we can create a (text) meta index file from those `tbi` files:

```
❯ ./multi_tbx create_index my_tbi_files > my_index
```

Then we can use `./multi-tbx dump_region` to get VCF records for a pre-specified region:
```
❯ ./multi_tbx dump_region my_index chr10:1,000,000-1,010,000 > output
```

The tool is useful for quick look up when the variant data is scattered in many files by collecting information in the tbi file in one place. It handles the simple but tedius-to-do-manually logic for fetch a set of variants in a region automatically.


Usage:

```
❯ multi_tbx --help
multi_tbx 0.1.0
Jason Chin


USAGE:
    multi_tbx [SUBCOMMAND]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    create_index
    dump_region
    help            Prints this message or the help of the given subcommand(s)
```

```
❯ multi_tbx create_index --help
multi_tbx-create_index

USAGE:
    multi_tbx create_index <tbi_files>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

ARGS:
    <tbi_files>    Path to a list of tabix index files
```

```
❯ multi_tbx dump_region --help
multi_tbx-dump_region

USAGE:
    multi_tbx dump_region [FLAGS] <index_file> <region>

FLAGS:
        --col                integer, optional, specific the column (default to the 2nd column) for coordinates
    -h, --help               Prints help information
        --only_file_path     just show the vcf.gz file locations
    -V, --version            Prints version information
        --use_whole_block    dump whole index block

ARGS:
    <index_file>    Path to a meta tabix index file
    <region>        the region of interest in the format {chr_str}:{bgn_u32}-{end_u32}
```
