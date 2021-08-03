# ==============================================================================
# Environment
# ==============================================================================
from interval import interval

# ==============================================================================
# Classes
# ==============================================================================
class GenomicInterval:
	def __init__(self, chrom, start, end, data=None):
		self.chr = chrom
		self.interval = interval([start, end])
		self.data = data

	def coord_str(self):
		return (
			self.chr
			+ " ["
			+ coord_pos_str(int(self.interval[0].inf))
			+ "e6, "
			+ coord_pos_str(int(self.interval[0].sup))
			+ "e6)"
		)

	def plot_str(self):
		return "".join(
			[
				self.chr,
				":",
				coord_pos_str(int(self.interval[0].inf)),
				"Mb-",
				coord_pos_str(int(self.interval[0].sup)),
				"Mb",
			]
		)

	def __str__(self):
		return self.plot_str()

	def __repr__(self):
		return self.plot_str()

	def inf(self):
		return int(self.interval[0].inf)

	def sup(self):
		return int(self.interval[0].sup)
	
	def __lt__(self, other):
		if self.chr == other.chr:
			return self.sup() < other.sup()
		return self.chr < other.chr


class TopologicalDomain(GenomicInterval):
	def __init__(self, chrom, start, end, persistence=(None, None)):
		self.persistence = persistence
		super().__init__(chrom, start, end)


# ==============================================================================
# Functions
# ==============================================================================


def coord_pos_str(i: int):
	return str(round(i / 1e6, 2))


def overlapping(a: GenomicInterval, b: GenomicInterval, extend: int = 0):
	"""
	Boolean function to see if two GenomicIntervals overlap

	Parameters
	----------
	a: GenomicInterval
	b: GenomicInterval
	extend: int
		Amount to extend each interval by before comparing
	"""
	return (
		a.chr == b.chr
		and a.inf() - extend <= b.sup() + extend
		and b.inf() - extend <= a.sup() + extend
	)


def get_mutated_ids_near_breakend(bp, neighbours, tol=1e5):
	"""
	Get the sample IDs of all samples with a breakpoint near a given breakpoint end

	Parameters
	----------
	bp : nx.Node
		Breakpoint under consideration
	neighbours : [nx.Node]
		Nearby and recurrent neighbours of `bp`
	tol : int
		Distance tolerance around `pos` to be considered for "recurrent"
	"""
	return list(
		set(
			[bp.data["sample"]]
			+ [n.data["sample"] for n in neighbours if overlapping(bp, n, tol)]
		)
	)


def find_tad(i: GenomicInterval, tads):
	"""
	Find the parent TAD(s) containing this genomic interval

	Parameters
	----------
	i : GenomicInterval
		Interval to find the TAD of
	tads : pd.DataFrame
		TADs for a given patient called at a given window size
	"""
	# get TAD containing infimum
	tads_lower = tads.loc[
		(tads.chr == i.chr) & (tads.start <= i.inf()) & (i.inf() <= tads.end), :
	]
	# get TAD containing supremum
	tads_upper = tads.loc[
		(tads.chr == i.chr) & (tads.start <= i.sup()) & (i.sup() <= tads.end), :
	]
	# get indices of the TADs identified
	lower_idx = tads_lower.index.tolist()
	upper_idx = tads_upper.index.tolist()
	# if everything is contained to a single TAD, return just the one
	if (
		(len(lower_idx) == 1)
		and (len(upper_idx) == 1)
		and (lower_idx[0] == upper_idx[0])
	):
		return tads_lower
	# if not, return the set of TADs spanned by the breakpoint coordinates
	else:
		return tads.iloc[range(min(lower_idx), max(upper_idx) + 1), :]

