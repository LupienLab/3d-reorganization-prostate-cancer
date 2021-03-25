# ==============================================================================
# Environment
# ==============================================================================
from interval import interval
from typing import Tuple

# ==============================================================================
# Classes
# ==============================================================================
class GenomicInterval:
    def __init__(self, chrom, start, end, data=None):
        self.chr = chrom
        self.interval = interval([start, end])
        self.data = data

    def coord_str(self) -> str:
        return (
            self.chr
            + " ["
            + coord_pos_str(int(self.interval[0].inf))
            + "e6, "
            + coord_pos_str(int(self.interval[0].sup))
            + "e6)"
        )

    def plot_str(self) -> str:
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

    def __str__(self) -> str:
        return self.plot_str()

    def __repr__(self) -> str:
        return self.plot_str()

    def inf(self) -> int:
        return int(self.interval[0].inf)

    def sup(self) -> int:
        return int(self.interval[0].sup)

    def __lt__(self, other) -> bool:
        if self.chr == other.chr:
            return self.inf() < other.inf()
        return self.chr < other.chr

    def __le__(self, other) -> bool:
        if self.chr == other.chr:
            return self.inf() <= other.inf()
        return self.chr <= other.chr

    def __gt__(self, other) -> bool:
        if self.chr == other.chr:
            return self.sup() > other.sup()
        return self.chr > other.chr

    def __ge__(self, other) -> bool:
        if self.chr == other.chr:
            return self.sup() >= other.sup()
        return self.chr >= other.chr

    def __eq__(self, other) -> bool:
        return (
            (self.chr == other.chr)
            and (self.inf() == other.inf())
            and (self.sup() == other.sup())
        )

    def __hash__(self):
        return hash(tuple((self.chr, self.inf(), self.sup(), tuple(self.data.items()))))


class TopologicalDomain(GenomicInterval):
    def __init__(self, chrom, start, end, persistence=(None, None)):
        self.persistence = persistence
        super().__init__(chrom, start, end)

    def __eq__(self, other) -> bool:
        return (
            (self.chr == other.chr)
            and (self.inf() == other.inf())
            and (self.sup() == other.sup())
            and (self.persistence == other.persistence)
        )


class Loop:
    def __init__(
        self,
        chrom_x,
        start_x,
        end_x,
        chrom_y,
        start_y,
        end_y,
        x_data=None,
        y_data=None,
        loop_data=None,
    ):
        # create anchors
        anchor_x = GenomicInterval(
            chrom_x,
            start_x,
            end_x,
            data=x_data,
        )
        anchor_y = GenomicInterval(
            chrom_y,
            start_y,
            end_y,
            data=y_data,
        )
        # anchor with the smaller coordinate is always the 0th element
        if anchor_x <= anchor_y:
            self.left = anchor_x
            self.right = anchor_y
        else:
            self.left = anchor_y
            self.right = anchor_x
        self.data = loop_data

    def gap(self) -> int:
        if self.left.chr != self.right.chr:
            return -1
        else:
            return self.right.inf() - self.left.sup()

    def __str__(self) -> str:
        return self.left.plot_str() + "--" + self.right.plot_str()

    def __repr__(self) -> str:
        return self.left.plot_str() + "--" + self.right.plot_str()

    def __lt__(self, other) -> bool:
        # use "less than" on the left anchor
        return self.left < other.left

    def __le__(self, other) -> bool:
        return self.left <= other.left

    def __gt__(self, other) -> bool:
        # use "greater than" on the right anchor
        return self.right > other.right

    def __ge__(self, other) -> bool:
        return self.right >= other.right

    def __eq__(self, other) -> bool:
        return (self.left == other.left) and (self.right == other.right)


# ==============================================================================
# Functions
# ==============================================================================


def coord_pos_str(i: int) -> str:
    return str(round(i / 1e6, 2))


def overlapping(a: GenomicInterval, b: GenomicInterval, extend: int = 0) -> bool:
    """
    Determine if two GenomicIntervals overlap

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
        range_slice = list(range(min(lower_idx), max(upper_idx) + 1))
        return GenomicInterval(
            i.chr,
            tads.loc[range_slice, "start"].min(),
            tads.loc[range_slice, "end"].max(),
        )
