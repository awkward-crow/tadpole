pub type ValueT = f32;
pub type IndexT = i64;
pub type CoefficientT = u16;

/// Format a float the same way C++ std::cout does for float: 6 significant digits,
/// trailing zeros stripped, scientific notation when exponent < -4 or >= 6.
pub fn fmt_float(v: ValueT) -> String {
    if v == 0.0 {
        return "0".to_string();
    }
    let exp = v.abs().log10().floor() as i32;
    if exp >= 6 || exp < -4 {
        // Scientific notation with 5 digits after decimal = 6 sig figs
        let s = format!("{:.5e}", v);
        // Rust produces e.g. "1.23450e6"; C++ produces "1.2345e+06"
        // Strip trailing zeros in the mantissa, then fix exponent formatting
        let (mantissa, exponent) = s.split_once('e').unwrap();
        let mantissa = mantissa.trim_end_matches('0').trim_end_matches('.');
        let exp_val: i32 = exponent.parse().unwrap();
        format!("{}e+{:02}", mantissa, exp_val)
    } else {
        // Fixed notation: enough decimal places for 6 significant digits
        let decimal_places = (5 - exp).max(0) as usize;
        let s = format!("{:.prec$}", v, prec = decimal_places);
        // Strip trailing zeros and trailing decimal point
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s
        }
    }
}

pub const NUM_COEFFICIENT_BITS: u32 = 8;

pub const MAX_SIMPLEX_INDEX: IndexT =
    (1_i64 << (8 * std::mem::size_of::<IndexT>() as u32 - 1 - NUM_COEFFICIENT_BITS)) - 1;

pub fn check_overflow(i: IndexT) {
    #[cfg(feature = "use_coefficients")]
    if i > MAX_SIMPLEX_INDEX {
        panic!(
            "simplex index {} in filtration is larger than maximum index {}",
            i as u64, MAX_SIMPLEX_INDEX
        );
    }
    #[cfg(not(feature = "use_coefficients"))]
    if i < 0 {
        panic!(
            "simplex index {} in filtration is larger than maximum index {}",
            i as u64, MAX_SIMPLEX_INDEX
        );
    }
}

// ---------------------------------------------------------------------------
// Entry — stores an index (and optionally a coefficient packed into the MSBs)
// ---------------------------------------------------------------------------

/// Without coefficients, Entry is just a wrapped IndexT.
/// With coefficients, the low NUM_COEFFICIENT_BITS bits hold the coefficient
/// and the remaining high bits hold the index.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Entry(pub IndexT);

impl Entry {
    #[cfg(not(feature = "use_coefficients"))]
    pub fn new(index: IndexT, _coefficient: CoefficientT) -> Self {
        Entry(index)
    }

    #[cfg(feature = "use_coefficients")]
    pub fn new(index: IndexT, coefficient: CoefficientT) -> Self {
        Entry((index << NUM_COEFFICIENT_BITS) | coefficient as IndexT)
    }

    #[cfg(not(feature = "use_coefficients"))]
    pub fn index(self) -> IndexT {
        self.0
    }

    #[cfg(feature = "use_coefficients")]
    pub fn index(self) -> IndexT {
        self.0 >> NUM_COEFFICIENT_BITS
    }

    #[cfg(not(feature = "use_coefficients"))]
    pub fn coefficient(self) -> CoefficientT {
        1
    }

    #[cfg(feature = "use_coefficients")]
    pub fn coefficient(self) -> CoefficientT {
        (self.0 & ((1 << NUM_COEFFICIENT_BITS) - 1)) as CoefficientT
    }

    #[cfg(feature = "use_coefficients")]
    pub fn set_coefficient(&mut self, c: CoefficientT) {
        self.0 = (self.0 & !((1 << NUM_COEFFICIENT_BITS) - 1)) | c as IndexT;
    }

    #[cfg(not(feature = "use_coefficients"))]
    pub fn set_coefficient(&mut self, _c: CoefficientT) {}
}

// ---------------------------------------------------------------------------
// DiameterIndex — (diameter, index) pair, used for simplex lists
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiameterIndex {
    pub diameter: ValueT,
    pub index: IndexT,
}

impl DiameterIndex {
    pub fn new(diameter: ValueT, index: IndexT) -> Self {
        DiameterIndex { diameter, index }
    }
}

// Ordering: larger diameter first, then smaller index first (matches C++)
impl PartialOrd for DiameterIndex {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for DiameterIndex {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other
            .diameter
            .partial_cmp(&self.diameter)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| self.index.cmp(&other.index))
    }
}

impl Eq for DiameterIndex {}

// ---------------------------------------------------------------------------
// IndexDiameter — (index, diameter) pair, used in sparse neighbor lists
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IndexDiameter {
    pub index: IndexT,
    pub diameter: ValueT,
}

impl IndexDiameter {
    pub fn new(index: IndexT, diameter: ValueT) -> Self {
        IndexDiameter { index, diameter }
    }
}

// Sorted by index for binary search in sparse distance matrix
impl PartialOrd for IndexDiameter {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for IndexDiameter {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.index.cmp(&other.index)
    }
}

impl Eq for IndexDiameter {}

// ---------------------------------------------------------------------------
// DiameterEntry — (diameter, Entry) pair, the main working unit
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
pub struct DiameterEntry {
    pub diameter: ValueT,
    pub entry: Entry,
}

impl DiameterEntry {
    pub fn new(diameter: ValueT, index: IndexT, coefficient: CoefficientT) -> Self {
        DiameterEntry {
            diameter,
            entry: Entry::new(index, coefficient),
        }
    }

    pub fn from_diameter_index(di: DiameterIndex, coefficient: CoefficientT) -> Self {
        DiameterEntry::new(di.diameter, di.index, coefficient)
    }

    pub fn invalid() -> Self {
        DiameterEntry::new(0.0, -1, 0)
    }

    pub fn index(self) -> IndexT {
        self.entry.index()
    }

    pub fn coefficient(self) -> CoefficientT {
        self.entry.coefficient()
    }

    pub fn set_coefficient(&mut self, c: CoefficientT) {
        self.entry.set_coefficient(c);
    }
}

impl PartialEq for DiameterEntry {
    fn eq(&self, other: &Self) -> bool {
        self.diameter == other.diameter && self.index() == other.index()
    }
}

impl Eq for DiameterEntry {}

// Priority queue ordering: larger diameter first, then smaller index first
impl PartialOrd for DiameterEntry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for DiameterEntry {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other
            .diameter
            .partial_cmp(&self.diameter)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| self.index().cmp(&other.index()))
    }
}
