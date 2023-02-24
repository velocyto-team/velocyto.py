from collections import defaultdict
from typing import Optional

from loguru import logger

from .constants import MATCH_INSIDE, MATCH_OVER3END, MATCH_OVER5END, MIN_FLANK
from .feature import Feature
from .read import Read
from .segment_match import SegmentMatch
from .transcript_model import TranscriptModel


class TransciptsIndex:
    __slots__ = ["transcipt_models", "tidx", "maxtidx"]
    """Search help class used to find the transcipt models that a read is spanning/contained into"""

    def __init__(self, trascript_models: list[TranscriptModel]) -> None:
        self.transcipt_models = trascript_models
        # self.transcipt_models.sort()
        self.tidx = 0  # index of the current interval
        self.maxtidx = len(trascript_models) - 1

    @property
    def scan_not_terminated(self) -> bool:
        """Return false when all the chromosome has been scanned"""
        return self.tidx < self.maxtidx

    def find_overlapping_trascript_models(self, read: Read) -> set[TranscriptModel]:
        """Finds all the Transcript models the Read overlaps with

        Args
        ----
        read: Read
            the read object to be analyzed

        Returns
        -------
        matched_transcripts: set of `TranscriptModel`
            TranscriptModel the read is overlapping with and values the kind of overlapping
            it is one of MATCH_INSIDE (1), MATCH_OVER5END (2), MATCH_OVER3END (4)

        """
        matched_transcripts: set[TranscriptModel] = set()
        if len(self.transcipt_models) == 0:
            logger.error(f"TransciptsIndex {self} contains no intervals")
            return matched_transcripts

        tmodel = self.transcipt_models[self.tidx]  # current transcript model
        # Move forward until we find the position we will never search left of again (because the reads are ordered)
        while self.scan_not_terminated and tmodel.ends_upstream_of(read):
            # move to the next interval
            self.tidx += 1
            tmodel = self.transcipt_models[self.tidx]

        # Loop trough the mapping segments of a read (e.g. just one of an internal exon, generally 2 for a splice or intron)
        for segment in read.segments:
            # Local search for each segment move a little forward (this just moves a coupple of intervals)
            i = self.tidx
            while i < self.maxtidx and tmodel.starts_upstream_of(segment):
                # matchtype = 0  # No match
                if tmodel.intersects(segment):
                    # NOTE: do we need to append all the model or the id would be enough?
                    matched_transcripts.add(tmodel)
                # move to the next transcript model
                i += 1
                tmodel = self.transcipt_models[i]
        return matched_transcripts


class FeatureIndex:
    """Search help class used to find the intervals that a read is spanning"""

    def __init__(self, ivls: Optional[list[Feature]] = None) -> None:
        self.ivls = ivls if ivls is not None else []
        self.ivls.sort()  # NOTE: maybe I am sorting twice check what I do upon creation
        self.iidx = 0  # index of the current interval
        self.maxiidx = len(ivls) - 1
        # NOTE needs to be changed to consider situations of identical intervals but different objects

    @property
    def last_interval_not_reached(self) -> bool:
        return self.iidx < self.maxiidx

    def reset(self) -> None:
        """It set the current feature to the first feature"""
        self.iidx = 0

    def has_ivls_enclosing(self, read: Read) -> bool:
        """Finds out if there are intervals that are fully containing all the read segments

        Args
        ----
        read: Read
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
            ivl = self.ivls[self.iidx]
            while i < self.maxiidx and ivl.doesnt_start_after(segment):
                matchtype = 0  # No match
                if ivl.contains(segment):
                    matchtype = MATCH_INSIDE
                if ivl.start_overlaps_with_part_of(
                    segment
                ):  # NOTE: should this be elif or it makes sense to allow both?
                    matchtype |= MATCH_OVER5END
                if ivl.end_overlaps_with_part_of(segment):  # NOTE: should this be elif or it makes sense to allow both?
                    matchtype |= MATCH_OVER3END

                segment_matchtype |= matchtype
                # move to the next interval
                i += 1
                ivl = self.ivls[i]

            # If one of the segments does not match inside a repeat we return false
            if segment_matchtype ^ MATCH_INSIDE:
                return False
        # If I arrive at this point of the code all the segments matched inside
        return True

    def mark_overlapping_ivls(self, read: Read) -> None:
        """Finds the overlap between Read and Features and mark intronic features if spanned

        Args
        ----
        read: Read
            the read object to be analyzed

        Returns
        -------
        Nothing, it marks the Feature object (is_validated = True) if there is evidence of exon-intron spanning
        """

        # ## TEMPORARY ######
        # dump_list = []
        ####################

        if len(self.ivls) == 0:
            return
        feature: Feature = self.ivls[self.iidx]  # current interval
        # Move forward until we find the position we will never search left of again (because the reads are ordered)
        while self.last_interval_not_reached and feature.ends_upstream_of(read):
            # move to the next interval
            self.iidx += 1
            feature = self.ivls[self.iidx]

        # Loop trough the mapping segments of a read (e.g. just one of an internal exon, generally 2 for a splice. for intron???)
        for segment in read.segments:
            # Local search for each segment move a little forward (this just moves a couple of intervals)
            i = self.iidx
            feature = self.ivls[self.iidx]
            while i < self.maxiidx and feature.doesnt_start_after(segment):
                # NOTE in this way the checks will be repeated for every clone of an interval in each transcript model
                # I might want to cache the check results of the previous feature
                if feature.kind == ord("i"):
                    if feature.end_overlaps_with_part_of(segment):
                        downstream_exon = feature.get_downstream_exon()
                        if downstream_exon.start_overlaps_with_part_of(segment):
                            feature.is_validated = True
                            # ## TEMPORARY ######
                            # dump_list.append((read, segment, feature, "r"))
                            ####################
                    if feature.start_overlaps_with_part_of(segment):
                        upstream_exon = feature.get_upstream_exon()
                        if upstream_exon.end_overlaps_with_part_of(segment):
                            feature.is_validated = True
                            # ## TEMPORARY ######
                            # dump_list.append((read, segment, feature, "l"))
                            ####################
                    # if feature.contains(segment):
                    #     pass  # here the intron is not validated !
                elif feature.kind != ord("e"):
                    raise ValueError(f"Unrecognized type of genomic feature {chr(feature.kind)}")

                #  move to the next interval
                i += 1
                feature = self.ivls[i]

        # ## TEMPORARY ######
        # return dump_list
        ####################

    def find_overlapping_ivls(self, read: Read) -> dict[TranscriptModel, list[SegmentMatch]]:
        """Finds the possible overlaps between Read and Features and return a 1 read derived mapping record

        Arguments
        ---------
        read: Read
            the read object to be analyzed

        Returns
        -------
        mapping_record: dict[TranscriptModel, list[SegmentMatch]]
            A record of the mappings by transcript model.
            Every entry contains a list of segment matches that in turn contains information on the segment and the feature

        Note
        ----
        - It is possible that a segment overalps at the same time an exon and an intron (spanning segment)
        - It is not possible that a segment overalps at the same time two exons. In that case the read is splitted
        into two segments and  the Read attribute `is_spliced == True`.
        - Notice that the name of the function might be confousing. if there is a non valid overallapping an empty mappign record will be return
        - Also notice that returning an empty mapping record will cause the suppression of the counting of the molecule

        """

        mapping_record: dict[TranscriptModel, list[SegmentMatch]] = defaultdict(list)

        if len(self.ivls) == 0:
            return mapping_record

        feature: Feature = self.ivls[self.iidx]  # current interval
        # Move forward until we find the position we will never search left of again (because the reads are ordered)
        while self.last_interval_not_reached and feature.ends_upstream_of(read):
            # move to the next interval
            self.iidx += 1
            feature = self.ivls[self.iidx]

        # Loop trough the mapping segments of a read (e.g. just one of an internal exon, generally 2 for a splice. for intron???)
        for segment in read.segments:
            # Local search for each segment move a little forward (this just moves a couple of intervals)
            i = self.iidx
            feature = self.ivls[i]
            while i < self.maxiidx and feature.doesnt_start_after(segment):
                # NOTE in this way the checks will be repeated for every clone of an interval in each transcript model
                # I might want to cache the check results of the previous feature
                # it was  if feature.contains(segment) or feature.start_overlaps_with_part_of(segment) or feature.end_overlaps_with_part_of(segment)
                # but I changed to the more simple
                if feature.intersects(segment) and (segment[-1] - segment[0]) > MIN_FLANK:
                    mapping_record[feature.transcript_model].append(SegmentMatch(segment, feature, read.is_spliced))
                #  move to the next interval
                i += 1
                feature = self.ivls[i]

        # NOTE: Removing first the one with less match and then requiring a splicing matching the transcript model is very stringent
        # It could be that for short ~10bp SKIP sequences the alligner has made a mistake and this might kill the whole molecule
        # NOTE: the code below is not very efficient
        if mapping_record:
            # Remove transcript models that are suboptimal match
            # in alternative one could use len(read.segment)
            max_n_segments = len(max(mapping_record.values(), key=len))
            for tm, segmatch_list in list(mapping_record.items()):
                if len(segmatch_list) < max_n_segments:
                    del mapping_record[tm]

        # NOTE: the code below is not very efficient, would be nice to avoid for loops
        # NOTE: potentailly bad effects: it could kill a lot of molecules if transcript models are not annotated correctly or missing
        if mapping_record:
            # A SKIP mapping needs to be explainaible by some kind of exon-exon exon-intron junction!
            # So if it falls internally, the TM needs to be removed from the mapping record
            for tm, segmatch_list in list(mapping_record.items()):
                for sm in segmatch_list:
                    if not sm.skip_makes_sense:
                        del mapping_record[tm]
                        break

        return mapping_record
