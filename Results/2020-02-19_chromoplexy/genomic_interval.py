# ==============================================================================
# Environment
# ==============================================================================
from interval import interval

# ==============================================================================
# Classes
# ==============================================================================
class GenomicInterval:
    def __init__(self, chr, start, end, data=None):
        self.chr = chr
        self.interval = interval([start, end])
        self.data = data

    def coord_str(self):
        return (
            self.chr
            + " ["
            + coord_pos_str(int(self.interval[0].inf))
            + ", "
            + coord_pos_str(int(self.interval[0].sup))
            + ")"
        )

    def plot_str(self):
        return "\n".join(
            [
                self.chr,
                coord_pos_str(int(self.interval[0].inf)),
                coord_pos_str(int(self.interval[0].sup)),
            ]
        )

    def __str__(self):
        return self.plot_str()

    def __repr__(self):
        return self.coord_str()

    def inf(self):
        return int(self.interval[0].inf)

    def sup(self):
        return int(self.interval[0].sup)


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
