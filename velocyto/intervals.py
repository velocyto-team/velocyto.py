from typing import *
from collections import defaultdict
import logging
import velocyto as vcy


class Interval:
    """ Holds an exon/intron interval read from the interval file """
    # Since many of this object are instanciated, implementing slots make it more efficient
    __slots__ = ['start', 'end', 'gene', "ivlidx", "ivltype", "is_sure_exon", "is_maybe_exon", "is_sure_intron", "is_sure_valid_intron", "is_last3prime"]

    def __init__(self, start: str, end: str, gene: vcy.Gene, ivlidx: int, ivltype: str, is_last3prime: bool) -> None:
        self.start, self.end, self.gene, self.ivlidx, self.ivltype = int(start), int(end), gene, ivlidx, ivltype
        self.is_last3prime = is_last3prime
        self.is_maybe_exon = ivltype in ('Eei', 'E-i')
        self.is_sure_exon = ivltype in ('E--', 'Ee-')
        self.is_sure_intron = ivltype in ('I--', 'I-i')
        self.is_sure_valid_intron = self.is_sure_intron  # NOTE This is just initialized and it will be adjusted later based on junction analysis
    
    def add_read_inside(self) -> None:
        self.gene.ivlinside_read_counts[self.ivlidx] += 1

    def add_read_spanning5end(self) -> None:
        self.gene.ivljunction5_read_counts[self.ivlidx] += 1

    def add_read_spanning3end(self) -> None:
        self.gene.ivljunction3_read_counts[self.ivlidx] += 1

    def ends_upstream_of(self, read: vcy.Read) -> bool:
        """The following situation happens
                                                            Read
                                               *|||segment|||-?-||segment|||????????
                ???????|||||Ivl|||||||||*

        """
        return self.end < read.pos

    def starts_upstream_of(self, segment: Tuple[int, int]) -> bool:
        """The following situation happens

                            *||||||segment|||||????????
                *|||||||||||||Ivl||||||||||????????????

        """
        return self.start < segment[-1]

    def contains(self, segment: Tuple[int]) -> bool:
        """The following situation happens

                     *||||||segment|||||*
                *|||||||||||||Ivl||||||||||||||||*

        """
        return (segment[0] >= self.start) and (segment[-1] <= self.end)

    def start_overlaps_with_part_of(self, segment: Tuple[int, int], minimum_flanking: int=vcy.MIN_FLANK) -> bool:
        """The following situation happens

          *---|||segment||---*
                *|||||||||||||Ivl||||||||||||||||*

        where `---` idicates the minimum flanking

        """
        return (segment[0] + minimum_flanking < self.start) and (segment[-1] - minimum_flanking > self.start)

    def end_overlaps_with_part_of(self, segment: Tuple[int, int], minimum_flanking: int=vcy.MIN_FLANK) -> bool:
        """The following situation happens

                                      *---|||segment||---*
                *|||||||||||||Ivl||||||||||||||||*

        where `---` idicates the minimum flanking
            
        """
        return (segment[0] + minimum_flanking < self.end) and (segment[-1] - minimum_flanking > self.end)

    def __lt__(self, other: Any) -> bool:
        assert type(other) == Interval, "Comparison can be only between two intervals"
        if self.start == other.start:
            return self.end < other.end
        return self.start < other.start

    def __gt__(self, other: Any) -> bool:
        assert type(other) == Interval, "Comparison can be only between two intervals"
        if self.start == other.start:
            return self.end > other.end
        return self.start > other.start


class IntervalsIndex:
    """ Search help class used to find the intervals that a read is spanning """
    def __init__(self, ivls: List[Interval]) -> None:
        self.ivls = ivls
        self.ivls.sort()  # NOTE: maybe I am sorting twice check what I do upon creation
        self.iidx = 0  # index of the current interval
        self.maxiidx = len(ivls) - 1

    @ property
    def last_interval_not_reached(self) -> bool:
        return self.iidx < self.maxiidx

    def has_ivls_enclosing(self, read: vcy.Read) -> bool:
        """Finds out if there are intervals that are fully containing all the read segments

        Args
        ----
        read: vcy.Read
            the read object to be analyzed

        Returns
        -------
        respones: bool
            if one has been found

        """
        if len(self.ivls) == 0:
            return False

        ivl = self.ivls[self.iidx]  # current interval
        
        # Move forward until we find the position we will never search left of again (because the reads are ordered)
        while self.last_interval_not_reached and ivl.ends_upstream_of(read):
            # move to the next interval
            self.iidx += 1
            ivl = self.ivls[self.iidx]
        
        for segment in read.segments:
            segment_matchtype = 0
            # Local search for each segment move a little forward (this just moves a coupple of intervals)
            i = self.iidx
            while i < self.maxiidx and ivl.starts_upstream_of(segment):
                matchtype = 0  # No match
                if ivl.contains(segment):
                    matchtype = vcy.MATCH_INSIDE
                if ivl.start_overlaps_with_part_of(segment):  # NOTE: should this be elif or it makes sense to allow both?
                    matchtype |= vcy.MATCH_OVER5END
                if ivl.end_overlaps_with_part_of(segment):  # NOTE: should this be elif or it makes sense to allow both?
                    matchtype |= vcy.MATCH_OVER3END
                
                segment_matchtype |= matchtype
                # move to the next interval
                i += 1
                ivl = self.ivls[i]
            
            # If one of the segments does not match inside a repeat we return false
            if segment_matchtype ^ vcy.MATCH_INSIDE:
                return False
        # If I arrive at this point of the code all the segment matched inside
        return True
        
    def find_overlapping_ivls(self, read: vcy.Read) -> Tuple[set, Dict[Any, int]]:
        """Finds the overlap between Read and intervals

        Args
        ----
        read: vcy.Read
            the read object to be analyzed

        Returns
        -------
        matchgenes: set
            the genes that the read is overlapping with

        matchivls: dict: {vcy.Interval: int}
            A dictionary with keys the intervals that the read is overlapping with and values the kind of overlapping
            it is one of vcy.MATCH_INSIDE (1), vcy.MATCH_OVER5END (2), vcy.MATCH_OVER3END (4)

        """
        matchgenes = set()  # type: set
        matchivls = defaultdict(int)  # type: Dict[Any, int]
        if len(self.ivls) == 0:
            # logging.warn("IntervalsIndex %s contains no intervals" % self)
            return matchgenes, matchivls

        ivl = self.ivls[self.iidx]  # current interval
        
        # Move forward until we find the position we will never search left of again (because the reads are ordered)
        while self.last_interval_not_reached and ivl.ends_upstream_of(read):
            # move to the next interval
            self.iidx += 1
            ivl = self.ivls[self.iidx]
        
        # Loop trough the mapping segments of a read (e.g. just one of an internal exon, generally 2 for a splice or intron)
        for segment in read.segments:
            
            # Local search for each segment move a little forward (this just moves a coupple of intervals)
            i = self.iidx
            while i < self.maxiidx and ivl.starts_upstream_of(segment):
                matchtype = 0  # No match
                if ivl.contains(segment):
                    matchtype = vcy.MATCH_INSIDE
                if ivl.start_overlaps_with_part_of(segment):  # NOTE: should this be elif or it makes sense to allow both?
                    matchtype |= vcy.MATCH_OVER5END
                if ivl.end_overlaps_with_part_of(segment):  # NOTE: should this be elif or it makes sense to allow both?
                    matchtype |= vcy.MATCH_OVER3END
                if matchtype:
                    # if there is a segment matching add it to the set of genes
                    matchgenes.add(ivl.gene)
                    # and add the match type to the set of genes
                    matchivls[ivl] |= matchtype
                # move to the next interval
                i += 1
                ivl = self.ivls[i]
    
        return matchgenes, matchivls
