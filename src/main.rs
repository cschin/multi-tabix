#[macro_use]
extern crate lazy_static;

use bgzip::tabix::{Tabix, TabixBin};
use bgzip::BGZFReader;
use clap::clap_app;
use regex::Regex;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug)]
struct IndexRecord {
    chr_name: String,
    bgn: usize,
    end: usize,
    bin: u32,
    file_path: String,
    chunk: u32,
    vfs_bgn: usize,
    vfs_end: usize,
}

fn dump_meta_index(path: String) -> Result<(), std::io::Error> {
    let f = BufReader::new(File::open(&path)?);
    f.lines()
        .into_iter()
        .try_for_each(|line| -> Result<(), std::io::Error> {
            let tbi_filename = line?;
            let tbi_filename = tbi_filename.trim_end_matches('\0');
            let mut tbi_file = File::open(tbi_filename)?;

            let tbx = Tabix::from_reader(&mut tbi_file)?;

            let idx = (0..37449)
                .into_iter()
                .map(|k| {
                    let l = ((7_f32 * k as f32 + 1 as f32).log2() / 3.0).floor() as u32;
                    let s = 2_u32.pow(29 - 3 * l) as u32;
                    let o = (8_u32.pow(l) - 1) / 7;
                    (k, ((k - o) * s, (k - o + 1) * s))
                })
                .collect::<FxHashMap<u32, (u32, u32)>>();

            (0..tbx.names.len()).into_iter().for_each(|i| {
                let seq_name = String::from_utf8_lossy(&tbx.names[i]);
                let seq_name = seq_name.trim_end_matches('\0');

                let mut r = tbx.sequences[i]
                    .bins
                    .iter()
                    .filter(|&x| *x.0 <= 37449)
                    .map(|x| (idx.get(x.0).unwrap().to_owned(), *x.0, x.1))
                    .collect::<Vec<((u32, u32), u32, &TabixBin)>>();
                r.sort_by_key(|k| k.0 .0);
                r.iter().for_each(|rr| {
                    rr.2.chunks.iter().enumerate().for_each(|(cid, chunk)| {
                        println!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            seq_name,
                            rr.0 .0,
                            rr.0 .1,
                            rr.1,
                            tbi_filename.strip_suffix(".tbi").unwrap(),
                            cid,
                            chunk.begin,
                            chunk.end
                        )
                    })
                });
            });
            Ok(())
        })?;
    Ok(())
}

fn read_meta_index(
    path: String,
) -> Result<FxHashMap<(String, u32), Vec<IndexRecord>>, std::io::Error> {
    let mut index = FxHashMap::<(String, u32), Vec<IndexRecord>>::default();
    let f = BufReader::new(File::open(&path)?);

    f.lines()
        .into_iter()
        .try_for_each(|line| -> Result<(), std::io::Error> {
            let line = line?;
            //println!("{}", line);
            let line = line.trim_end_matches("\n\r");
            let fields: Vec<&str> = line.split('\t').collect();
            let rec = IndexRecord {
                chr_name: fields[0].to_owned(),
                bgn: fields[1]
                    .parse::<usize>()
                    .expect(format!("parsing error: bgn = {}", fields[1]).as_str()),
                end: fields[2]
                    .parse::<usize>()
                    .expect(format!("parsing error: end = {}", fields[2]).as_str()),
                bin: fields[3]
                    .parse::<u32>()
                    .expect(format!("parsing error: bin = {}", fields[3]).as_str()),
                file_path: String::from(fields[4]),
                chunk: fields[5]
                    .parse::<u32>()
                    .expect(format!("parsing error: chunk = {}", fields[5]).as_str()),
                vfs_bgn: fields[6]
                    .parse::<usize>()
                    .expect(format!("parsing error: vfs_bgn = {}", fields[6]).as_str()),
                vfs_end: fields[7]
                    .parse::<usize>()
                    .expect(format!("parsing error: vfs_end = {}", fields[7]).as_str()),
            };
            let chr_name = fields[0].to_string();
            let bin = rec.bin;
            let key = (chr_name, bin);
            let e = index.entry(key).or_insert_with(|| vec![]);
            e.push(rec);
            Ok(())
        })?;
    Ok(index)
}

/*
static inline int reg2bins(uint32_t beg, uint32_t end, uint16_t list[MAX_BIN])
{
    int i = 0, k;
    if (beg >= end) return 0;
    if (end >= 1u<<29) end = 1u<<29;
    --end;
    list[i++] = 0;
    for (k =    1 + (beg>>26); k <=    1 + (end>>26); ++k) list[i++] = k;
    for (k =    9 + (beg>>23); k <=    9 + (end>>23); ++k) list[i++] = k;
    for (k =   73 + (beg>>20); k <=   73 + (end>>20); ++k) list[i++] = k;
    for (k =  585 + (beg>>17); k <=  585 + (end>>17); ++k) list[i++] = k;
    for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list[i++] = k;
    return i;
}
*/

fn rng2bins(bgn: u32, end: u32) -> Vec<u32> {
    let mut bins = Vec::<u32>::new();
    assert!(bgn < end);
    assert!(end < 1_u32 << 29);
    let mut end = end;
    end -= 1;
    bins.push(0);
    for k in (1 + (bgn >> 26))..=(1 + (end >> 26)) {
        bins.push(k);
    }
    for k in (9 + (bgn >> 23))..=(9 + (end >> 23)) {
        bins.push(k);
    }
    for k in (73 + (bgn >> 20))..=(73 + (end >> 20)) {
        bins.push(k);
    }
    for k in (585 + (bgn >> 17))..=(585 + (end >> 17)) {
        bins.push(k);
    }
    for k in (4681 + (bgn >> 14))..=(4681 + (end >> 14)) {
        bins.push(k);
    }
    bins
}

// TODO: parsing error handling
fn parse_rgn(s: String) -> (String, u32, u32) {
    lazy_static! {
        static ref RGN: Regex = Regex::new(r"^([a-zA-Z0-9]+):(\d+)-(\d+)$").unwrap();
    }
    let s = s.replace(" ", "").replace(",", "");
    let caps = RGN
        .captures(&s)
        .expect(format!("region format error: {}", s).as_str());
    //println!("{:?}", caps);
    let seq_name = caps
        .get(1)
        .expect(format!("seq name error").as_str())
        .as_str()
        .to_string();
    let bgn = caps
        .get(2)
        .expect(format!("can'get the reg-begin").as_str())
        .as_str()
        .parse::<u32>()
        .expect("reg-begin parsing error");
    let end = caps
        .get(3)
        .expect(format!("can'get the reg-end").as_str())
        .as_str()
        .parse::<u32>()
        .expect("reg-end parsing error");

    (seq_name, bgn, end)
}

fn main() -> Result<(), std::io::Error> {
    let app_m = clap_app!(multi_tbx =>
        (version: "0.1.0")
        (author: "Jason Chin")
        (about: "")
        (@subcommand create_index =>
            (@arg tbi_files: +required "Path to a list of tabix index files")
        )
        (@subcommand dump_region =>
            (@arg index_file: +required "Path to a meta tabix index file")
            (@arg region: +required "the region of interest in the format {chr_str}:{bgn_u32}-{end_u32}")
            (@arg whole_block: --use_whole_block "dump whole index block")
            (@arg only_file_path: --only_file_path "just show the vcf.gz file locations")
            (@arg coordinate_column: --col "integer, optional, specific the column (default to the 2nd column) for coordinates")
        )
    )
    .get_matches();

    match app_m.subcommand() {
        ("create_index", Some(sub_m)) => {
            let tbi_files = sub_m
                .value_of("tbi_files")
                .expect("error to get the list of tbi files")
                .to_string();
            dump_meta_index(tbi_files)?;
        }
        ("dump_region", Some(sub_m)) => {
            let index_file = sub_m
                .value_of("index_file")
                .expect("error to get the index_file name")
                .to_string();
            let rgn_spec = sub_m
                .value_of("region")
                .expect("error to get the region")
                .to_string();

            let use_whole_block = sub_m.is_present("whole_block");
            let only_file_path = sub_m.is_present("only_file_path");
            let coordinate_col = match sub_m.value_of("coordinate_column") {
                Some(col) => {
                    col.parse::<u32>()
                        .expect(format!("parsing error for col {}", col).as_str())
                        - 1
                }
                _ => 1,
            };

            let index_rec = read_meta_index(index_file)?;
            let rgn = parse_rgn(rgn_spec);
            let target_bins = rng2bins(rgn.1, rgn.2);

            let mut vfs_offsets = FxHashMap::<String, Vec<(usize, usize)>>::default();
            target_bins.iter().for_each(|&bin| {
                //println!("bin: {}", bin);
                let seq_name = rgn.0.to_owned();
                let key = (seq_name, bin);
                if index_rec.contains_key(&key) {
                    let recs = index_rec.get(&key).unwrap();
                    recs.iter().for_each(|r| {
                        let e = vfs_offsets
                            .entry(r.file_path.to_owned())
                            .or_insert_with(|| vec![]);
                        e.push((r.vfs_bgn, r.vfs_end));
                    });
                }
            });

            vfs_offsets.iter_mut().for_each(|(path, offsets)| {

                if only_file_path {
                    println!("{}", path);
                    return;
                }

                offsets.sort();
                let mut reader = BGZFReader::new(File::open(&path).unwrap());
                let mut line = String::new();
                let cur_offset = offsets[0].0 as u64;
                let max_offset = offsets[offsets.len() - 1].1 as u64;
                reader.bgzf_seek(cur_offset as u64).unwrap();
           

                if use_whole_block {
                    loop {
                        line.clear();
                        reader.read_line(&mut line).expect("error on getting vcf record");
                        print!("{}", line);
                        if reader.bgzf_pos() > max_offset {
                            break;
                        }
                    }
                } else {
                    let mut cur_pos = 0_u32;
                    loop {
                        line.clear();
                        reader.read_line(&mut line).expect("error on getting vcf record");
                        let splitted = line.split_whitespace();
                        for (c_idx, s) in splitted.enumerate() {
                            if coordinate_col == c_idx as u32 {
                                let coordinate = s
                                    .parse::<u32>()
                                    .expect(format!("coordinate parsing error: {}", s).as_str());
                                cur_pos = coordinate;
                                if coordinate >= rgn.1 && coordinate <= rgn.2 {
                                    print!("{}", line);
                                }
                                break;
                            }
                            if coordinate_col < c_idx as u32 {
                                continue;
                            };
                        };

                        if cur_pos > rgn.2 {
                            break;
                        }
                        if reader.bgzf_pos() > max_offset {
                            break;
                        }
                    }
                }
            });
        }
        _ => {}
    };
    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_parse_rgn() {
        let out = crate::parse_rgn("test:1000-50000".to_string());
        assert_eq!(out, ("test".to_string(), 1000_u32, 50000_u32));
        let out = crate::parse_rgn("test :1,000 - 50,000".to_string());
        assert_eq!(out, ("test".to_string(), 1000_u32, 50000_u32));
    }
}
