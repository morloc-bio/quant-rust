// Rust implementations of quant, over salmon 2.7.0's eqclass, model and
// infer crates.

#[derive(Clone, Debug, Default)]
pub struct BqFragLenDist {
    pub counts: Vec<i64>,
}

#[derive(Clone, Debug, Default)]
pub struct BqStrandCounts {
    pub sense: i64,
    pub antisense: i64,
    pub mixed: i64,
}

#[allow(non_snake_case)]
#[derive(Clone, Debug)]
pub struct BqQuantModel {
    pub layout: BioLayout,
    pub strand: BioLibStrand,
    pub fldLogPmf: Vec<f64>,
    pub effLens: Vec<f64>,
}

#[derive(Clone, Debug, Default)]
pub struct BqEqTable {
    pub txs: Vec<Vec<i64>>,
    pub bins: Vec<Vec<i64>>,
    pub weights: Vec<Vec<f64>>,
    pub counts: Vec<i64>,
}

#[allow(non_snake_case)]
#[derive(Clone, Debug)]
pub struct BqEmOpts {
    pub vbem: bool,
    pub vbPrior: f64,
    pub maxIter: i64,
}

#[allow(non_snake_case)]
#[derive(Clone, Debug, Default)]
pub struct BqTxQuant {
    pub counts: Vec<f64>,
    pub tpm: Vec<f64>,
    pub effLens: Vec<f64>,
}

#[derive(Clone, Debug, Default)]
pub struct BqGeneQuant {
    pub counts: Vec<f64>,
    pub tpm: Vec<f64>,
    pub lengths: Vec<f64>,
}

#[repr(u8)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum BqCountScaling { NoScaling = 0, ScaledTpm = 1, LengthScaledTpm = 2 }

// salmon's defaults.
const BQ_FLD_MAX: usize = 1000;
const BQ_FLD_MEAN: f64 = 250.0;
const BQ_FLD_SD: f64 = 25.0;
const BQ_RANGE_BINS: u32 = 4;

fn bq_add_into(a: &mut Vec<i64>, b: &[i64]) {
    if a.len() < b.len() {
        a.resize(b.len(), 0);
    }
    for (x, y) in a.iter_mut().zip(b) {
        *x += y;
    }
}

pub fn bq_frag_len_dist(_anchor: i64, xs: &Vec<BioFragHits>) -> BqFragLenDist {
    let mut counts = vec![0i64; BQ_FLD_MAX + 1];
    for f in xs {
        let Some(h0) = f.hits.first() else { continue };
        if h0.mates != 0 || h0.fragLen <= 0 {
            continue;
        }
        if f.hits.iter().all(|h| h.mates == 0 && h.fragLen == h0.fragLen) {
            counts[(h0.fragLen as usize).min(BQ_FLD_MAX)] += 1;
        }
    }
    BqFragLenDist { counts }
}

pub fn bq_merge_frag_len_dist(a: &BqFragLenDist, b: &BqFragLenDist) -> BqFragLenDist {
    let mut counts = a.counts.clone();
    bq_add_into(&mut counts, &b.counts);
    BqFragLenDist { counts }
}

// Whether read 1 lies on the transcript's forward strand: an orphaned mate 2
// reports its own strand.
fn bq_read1_fwd(h: &BioTxHit) -> bool {
    if h.mates == 2 { !h.fwd } else { h.fwd }
}

pub fn bq_strand_counts(_anchor: i64, xs: &Vec<BioFragHits>) -> BqStrandCounts {
    let mut c = BqStrandCounts::default();
    for f in xs {
        let Some(h0) = f.hits.first() else { continue };
        let s = bq_read1_fwd(h0);
        if f.hits.iter().all(|h| bq_read1_fwd(h) == s) {
            if s { c.sense += 1 } else { c.antisense += 1 }
        } else {
            c.mixed += 1;
        }
    }
    c
}

pub fn bq_merge_strand_counts(a: &BqStrandCounts, b: &BqStrandCounts) -> BqStrandCounts {
    BqStrandCounts { sense: a.sense + b.sense, antisense: a.antisense + b.antisense, mixed: a.mixed + b.mixed }
}

pub fn bq_strand_call(stranded: f64, unstranded: f64, c: &BqStrandCounts) -> String {
    let total = (c.sense + c.antisense) as f64;
    if total <= 0.0 {
        return "undetermined".to_string();
    }
    let (sf, sr) = (c.sense as f64 / total, c.antisense as f64 / total);
    let call = if sf >= stranded {
        "forward"
    } else if sr >= stranded {
        "reverse"
    } else if (sf - sr).abs() <= unstranded {
        "unstranded"
    } else {
        "undetermined"
    };
    call.to_string()
}

pub fn bq_quant_model(_anchor: i64, layout: BioLayout, strand: BioLibStrand, d: &BqFragLenDist, ref_lens: &Vec<i64>) -> BqQuantModel {
    let mut fld = salmon_model::FragmentLengthDistribution::new(1.0, BQ_FLD_MAX, BQ_FLD_MEAN, BQ_FLD_SD, 4, 0.5, 1);
    for (len, &c) in d.counts.iter().enumerate() {
        if c > 0 {
            fld.add_val(len, (c as f64).ln());
        }
    }
    fld.cache();
    let pmf: Vec<f64> = (0..=BQ_FLD_MAX).map(|i| fld.pmf(i)).collect();
    let cm = fld.conditional_means();
    let eff = ref_lens.iter().map(|&l| salmon_model::smoothed_effective_length(&cm, l.max(0) as usize)).collect();
    BqQuantModel { layout, strand, fldLogPmf: pmf, effLens: eff }
}

fn bq_expected(l: BioLayout, s: BioLibStrand) -> salmon_core::LibraryFormat {
    let code = match (l, s) {
        (BioLayout::PairedEnd, BioLibStrand::LibUnstranded) => "IU",
        (BioLayout::PairedEnd, BioLibStrand::LibForward) => "ISF",
        (BioLayout::PairedEnd, BioLibStrand::LibReverse) => "ISR",
        (BioLayout::SingleEnd, BioLibStrand::LibUnstranded) => "U",
        (BioLayout::SingleEnd, BioLibStrand::LibForward) => "SF",
        (BioLayout::SingleEnd, BioLibStrand::LibReverse) => "SR",
    };
    salmon_core::LibraryFormat::parse(code).unwrap_or_else(|e| panic!("library format {}: {:?}", code, e))
}

fn bq_status(l: BioLayout, mates: i64) -> salmon_core::MateStatus {
    match (l, mates) {
        (BioLayout::SingleEnd, _) => salmon_core::MateStatus::SingleEnd,
        (_, 1) => salmon_core::MateStatus::PairedEndLeft,
        (_, 2) => salmon_core::MateStatus::PairedEndRight,
        _ => salmon_core::MateStatus::PairedEndPaired,
    }
}

type BqClassMap = std::collections::HashMap<(Vec<u32>, Vec<u32>), (Vec<f64>, i64)>;

fn bq_eq_chunk(m: &BqQuantModel, xs: &[BioFragHits]) -> BqClassMap {
    let expected = bq_expected(m.layout, m.strand);
    let mut map = BqClassMap::new();
    let mut aux: Vec<(u32, f64)> = Vec::new();
    for f in xs {
        aux.clear();
        for h in &f.hits {
            let observed = if h.format >= 0 { Some(salmon_core::LibraryFormat::from_format_id(h.format as u8)) } else { None };
            if !salmon_core::is_compatible(expected, observed, h.fwd, bq_status(m.layout, h.mates)) {
                continue;
            }
            let mut lp = h.weight.max(f64::MIN_POSITIVE).ln();
            if h.mates == 0 && h.fragLen > 0 {
                lp += m.fldLogPmf[(h.fragLen as usize).min(m.fldLogPmf.len() - 1)];
            }
            aux.push((h.tx as u32, lp));
        }
        if aux.is_empty() {
            continue;
        }
        let mx = aux.iter().map(|x| x.1).fold(f64::NEG_INFINITY, f64::max);
        aux.sort_by_key(|x| x.0);
        let mut txs: Vec<u32> = Vec::with_capacity(aux.len());
        let mut ws: Vec<f64> = Vec::with_capacity(aux.len());
        for &(t, lp) in aux.iter() {
            let w = (lp - mx).exp();
            if txs.last() == Some(&t) {
                *ws.last_mut().unwrap() += w;
            } else {
                txs.push(t);
                ws.push(w);
            }
        }
        let tot: f64 = ws.iter().sum();
        for w in ws.iter_mut() {
            *w /= tot;
        }
        let bins = salmon_eqclass::range_factorize_bins(&ws, BQ_RANGE_BINS);
        let e = map.entry((txs, bins)).or_insert_with(|| (vec![0.0; ws.len()], 0));
        for (a, w) in e.0.iter_mut().zip(&ws) {
            *a += w;
        }
        e.1 += 1;
    }
    map
}

fn bq_merge_maps(mut a: BqClassMap, b: BqClassMap) -> BqClassMap {
    for (k, (ws, c)) in b {
        let e = a.entry(k).or_insert_with(|| (vec![0.0; ws.len()], 0));
        for (x, w) in e.0.iter_mut().zip(&ws) {
            *x += w;
        }
        e.1 += c;
    }
    a
}

fn bq_table(map: BqClassMap) -> BqEqTable {
    let mut rows: Vec<((Vec<u32>, Vec<u32>), (Vec<f64>, i64))> = map.into_iter().collect();
    rows.sort_by(|a, b| a.0.cmp(&b.0));
    let mut t = BqEqTable::default();
    for ((txs, bins), (ws, c)) in rows {
        t.txs.push(txs.into_iter().map(|x| x as i64).collect());
        t.bins.push(bins.into_iter().map(|x| x as i64).collect());
        t.weights.push(ws);
        t.counts.push(c);
    }
    t
}

fn bq_map(t: &BqEqTable) -> BqClassMap {
    let mut m = BqClassMap::with_capacity(t.counts.len());
    for i in 0..t.counts.len() {
        let k = (t.txs[i].iter().map(|&x| x as u32).collect(), t.bins[i].iter().map(|&x| x as u32).collect());
        m.insert(k, (t.weights[i].clone(), t.counts[i]));
    }
    m
}

pub fn bq_eq_table(_anchor: i64, m: &BqQuantModel, xs: &Vec<BioFragHits>) -> BqEqTable {
    bq_table(bq_eq_chunk(m, xs))
}

pub fn bq_merge_eq_table(a: &BqEqTable, b: &BqEqTable) -> BqEqTable {
    bq_table(bq_merge_maps(bq_map(a), bq_map(b)))
}

pub fn bq_quantify(o: &BqEmOpts, m: &BqQuantModel, t: &BqEqTable) -> BqTxQuant {
    let n = m.effLens.len();
    let builder = salmon_eqclass::EquivalenceClassBuilder::new();
    for i in 0..t.counts.len() {
        let txs: Vec<u32> = t.txs[i].iter().map(|&x| x as u32).collect();
        let bins: Vec<u32> = t.bins[i].iter().map(|&x| x as u32).collect();
        if txs.iter().any(|&x| x as usize >= n) {
            panic!("quantify: a class names reference {} of {}", txs.iter().max().unwrap(), n);
        }
        builder.add_group(salmon_eqclass::TranscriptGroup::with_bins(txs, bins), t.weights[i].clone(), t.counts[i] as u64);
    }
    let mut eq = builder.finish();
    eq.update_eff_lengths(&m.effLens);
    let mut opts = salmon_infer::EmOptions::default();
    opts.use_vbem = o.vbem;
    opts.vb_prior = o.vbPrior;
    opts.max_iter = o.maxIter.max(1) as u32;
    let r = salmon_infer::optimize(&eq, n, &opts, salmon_infer::EffLens::new(&m.effLens));
    let counts = r.alphas;
    let rates: Vec<f64> = counts.iter().zip(&m.effLens).map(|(c, l)| c / l.max(1.0)).collect();
    let tot: f64 = rates.iter().sum();
    let tpm = rates.iter().map(|r| if tot > 0.0 { 1e6 * r / tot } else { 0.0 }).collect();
    BqTxQuant { counts, tpm, effLens: m.effLens.clone() }
}

pub fn bq_gene_quant(gene_of: &Vec<i64>, ngenes: i64, q: &BqTxQuant) -> BqGeneQuant {
    let n = ngenes.max(0) as usize;
    let mut g = BqGeneQuant { counts: vec![0.0; n], tpm: vec![0.0; n], lengths: vec![0.0; n] };
    let mut lsum = vec![0.0; n];
    let mut ntx = vec![0usize; n];
    for (i, &gi) in gene_of.iter().enumerate() {
        if gi < 0 {
            continue;
        }
        let gi = gi as usize;
        g.counts[gi] += q.counts[i];
        g.tpm[gi] += q.tpm[i];
        g.lengths[gi] += q.tpm[i] * q.effLens[i];
        lsum[gi] += q.effLens[i];
        ntx[gi] += 1;
    }
    for i in 0..n {
        g.lengths[i] = if g.tpm[i] > 0.0 {
            g.lengths[i] / g.tpm[i]
        } else if ntx[i] > 0 {
            lsum[i] / ntx[i] as f64
        } else {
            0.0
        };
    }
    g
}

pub fn bq_scaled_counts(s: BqCountScaling, qs: &Vec<BqGeneQuant>) -> Vec<Vec<f64>> {
    if s == BqCountScaling::NoScaling || qs.is_empty() {
        return qs.iter().map(|q| q.counts.clone()).collect();
    }
    let n = qs[0].tpm.len();
    let mean_len: Vec<f64> = (0..n).map(|g| qs.iter().map(|q| q.lengths[g]).sum::<f64>() / qs.len() as f64).collect();
    qs.iter()
        .map(|q| {
            let raw: Vec<f64> = if s == BqCountScaling::LengthScaledTpm {
                q.tpm.iter().zip(&mean_len).map(|(t, l)| t * l).collect()
            } else {
                q.tpm.clone()
            };
            let want: f64 = q.counts.iter().sum();
            let have: f64 = raw.iter().sum();
            raw.iter().map(|x| if have > 0.0 { x * want / have } else { 0.0 }).collect()
        })
        .collect()
}

pub fn bq_format_id(code: &String) -> i64 {
    salmon_core::LibraryFormat::parse(code).unwrap_or_else(|e| panic!("library format {}: {:?}", code, e)).format_id() as i64
}

// salmon-shaped reports
// ---------------------

fn bq_fld_moments(fld: &BqFragLenDist) -> (f64, f64) {
    let n: i64 = fld.counts.iter().sum();
    if n == 0 {
        return (BQ_FLD_MEAN, BQ_FLD_SD);
    }
    let mean = fld.counts.iter().enumerate().map(|(i, &c)| i as f64 * c as f64).sum::<f64>() / n as f64;
    let var = fld.counts.iter().enumerate().map(|(i, &c)| c as f64 * (i as f64 - mean).powi(2)).sum::<f64>() / n as f64;
    (mean, var.sqrt())
}

fn bq_json_str(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

pub fn bq_salmon_meta_info(lib: &String, fragments: i64, mapped: i64, decoy: i64, orphan: i64, fld: &BqFragLenDist) -> String {
    let (mean, sd) = bq_fld_moments(fld);
    let pct = if fragments > 0 { 100.0 * mapped as f64 / fragments as f64 } else { 0.0 };
    format!(
        "{{\n  \"salmon_version\": \"biolib\",\n  \"samp_type\": \"none\",\n  \"opt_type\": \"vb\",\n  \"num_libraries\": 1,\n  \"library_types\": [{lib}],\n  \"frag_dist_length\": {len},\n  \"frag_length_mean\": {mean},\n  \"frag_length_sd\": {sd},\n  \"frag_length_source\": \"reads\",\n  \"seq_bias_correct\": false,\n  \"gc_bias_correct\": false,\n  \"pos_bias_correct\": false,\n  \"mapping_type\": \"mapping\",\n  \"num_processed\": {fragments},\n  \"num_mapped\": {mapped},\n  \"percent_mapped\": {pct},\n  \"num_decoy_fragments\": {decoy},\n  \"num_orphan\": {orphan},\n  \"num_bootstraps\": 0,\n  \"detected_library_type\": {lib},\n  \"call\": \"quant\"\n}}\n",
        lib = bq_json_str(lib),
        len = BQ_FLD_MAX + 1,
    )
}

pub fn bq_lib_format_counts(lib: &String, s: &BqStrandCounts) -> String {
    let both = s.sense + s.antisense;
    let compatible = match lib.as_str() {
        "ISF" | "SF" => s.sense,
        "ISR" | "SR" => s.antisense,
        "IU" | "U" => both,
        other => panic!("libFormatCounts: unknown library type '{}'", other),
    };
    let pre = if lib.starts_with('I') { "I" } else { "" };
    let ratio = if both > 0 { compatible as f64 / both as f64 } else { 0.0 };
    let bias = if both > 0 { s.sense as f64 / both as f64 } else { 0.0 };
    format!(
        "{{\n  \"expected_format\": {},\n  \"compatible_fragment_ratio\": {},\n  \"num_compatible_fragments\": {},\n  \"num_incompatible_fragments\": {},\n  \"num_assigned_fragments\": {},\n  \"num_frags_with_inconsistent_or_orphan_mappings\": {},\n  \"strand_mapping_bias\": {},\n  \"{pre}U\": 0,\n  \"{pre}SF\": {},\n  \"{pre}SR\": {}\n}}\n",
        bq_json_str(lib), ratio, compatible, both - compatible, both, s.mixed, bias, s.sense, s.antisense
    )
}

pub fn bq_flen_dist(fld: &BqFragLenDist) -> String {
    let n: i64 = fld.counts.iter().sum();
    let mut v = vec![0f64; BQ_FLD_MAX + 1];
    for (i, &c) in fld.counts.iter().enumerate() {
        if n > 0 {
            v[i.min(BQ_FLD_MAX)] += c as f64 / n as f64;
        }
    }
    let mut out = v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("\t");
    out.push('\n');
    out
}

pub fn bq_infer_experiment_text(layout: BioLayout, s: &BqStrandCounts) -> String {
    let n = (s.sense + s.antisense + s.mixed).max(1) as f64;
    let (kind, sense, anti) = match layout {
        BioLayout::PairedEnd => ("PairEnd", "1++,1--,2+-,2-+", "1+-,1-+,2++,2--"),
        BioLayout::SingleEnd => ("SingleEnd", "++,--", "+-,-+"),
    };
    format!(
        "\n\nThis is {} Data\nFraction of reads failed to determine: {:.4}\nFraction of reads explained by \"{}\": {:.4}\nFraction of reads explained by \"{}\": {:.4}\n",
        kind, s.mixed as f64 / n, sense, s.sense as f64 / n, anti, s.antisense as f64 / n
    )
}
