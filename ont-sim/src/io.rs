use crate::{Simulation, Truth};
use flate2::{Compression, write::GzEncoder};
use std::{
    fs::File,
    io::{self, BufRead, BufWriter, Write},
    path::Path,
    str::FromStr,
};

pub const TRUTH_HEADER: &str = "read_idx\tread_name\tmolecule_id\tgenome_id\tmol_start\tmol_end\tis_rc\tfoldback_len\tpolya_len\tn_sub\tn_ins\tn_del\tmean_q";

pub fn write_outputs(sim: &Simulation, fastq: &Path, truth: &Path) -> io::Result<()> {
    if fastq == truth {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "FASTQ and truth paths must differ",
        ));
    }
    let mut truth = BufWriter::new(File::create(truth)?);
    writeln!(truth, "{TRUTH_HEADER}")?;
    let mut file = BufWriter::new(File::create(fastq)?);
    let mut emit = |out: &mut dyn Write| -> io::Result<()> {
        for r in &sim.reads {
            let t = &r.truth;
            writeln!(out, "@{}", t.read_name)?;
            out.write_all(&r.sequence)?;
            out.write_all(b"\n+\n")?;
            out.write_all(&r.quality)?;
            out.write_all(b"\n")?;
            writeln!(
                truth,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                t.read_idx,
                t.read_name,
                t.molecule_id,
                t.genome_id,
                t.mol_start,
                t.mol_end,
                t.is_rc,
                t.foldback_len,
                t.polya_len,
                t.n_sub,
                t.n_ins,
                t.n_del,
                t.mean_q
            )?;
        }
        out.flush()?;
        truth.flush()
    };
    if fastq.to_string_lossy().ends_with(".gz") {
        let mut gz = GzEncoder::new(file, Compression::default());
        emit(&mut gz)?;
        gz.finish()?.flush()
    } else {
        emit(&mut file)
    }
}

fn parse<T: FromStr>(s: &str) -> Result<T, String> {
    s.parse().map_err(|_| format!("invalid value: {s:?}"))
}
fn table<T>(
    reader: impl BufRead,
    header: &str,
    parse: impl Fn(&[&str]) -> Result<T, String>,
) -> Result<Vec<T>, String> {
    let mut lines = reader.lines();
    if lines
        .next()
        .transpose()
        .map_err(|e| e.to_string())?
        .as_deref()
        != Some(header)
    {
        return Err(format!("expected TSV header: {header}"));
    }
    let columns = header.split('\t').count();
    lines
        .enumerate()
        .map(|(i, line)| {
            let line = line.map_err(|e| e.to_string())?;
            let fields: Vec<_> = line.split('\t').collect();
            let row = if fields.len() == columns {
                parse(&fields)
            } else {
                Err("wrong column count".into())
            };
            row.map_err(|e| format!("line {}: {e}", i + 2))
        })
        .collect()
}

pub fn read_truth(reader: impl BufRead) -> Result<Vec<Truth>, String> {
    table(reader, TRUTH_HEADER, |f| {
        let mean_q: f64 = parse(f[12])?;
        if f[1].is_empty() || !mean_q.is_finite() || !(0.0..=50.0).contains(&mean_q) {
            return Err("invalid read name or mean quality".into());
        }
        Ok(Truth {
            read_idx: parse(f[0])?,
            read_name: f[1].into(),
            molecule_id: parse(f[2])?,
            genome_id: parse(f[3])?,
            mol_start: parse(f[4])?,
            mol_end: parse(f[5])?,
            is_rc: parse(f[6])?,
            foldback_len: parse(f[7])?,
            polya_len: parse(f[8])?,
            n_sub: parse(f[9])?,
            n_ins: parse(f[10])?,
            n_del: parse(f[11])?,
            mean_q,
        })
    })
}

pub fn read_assignments(reader: impl BufRead) -> Result<Vec<(usize, usize)>, String> {
    table(reader, "read_idx\texemplar_idx", |f| {
        Ok((parse(f[0])?, parse(f[1])?))
    })
}
