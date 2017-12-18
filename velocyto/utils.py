from typing import *
import velocyto as vcy


def jump_next_3p_exon(feature: vcy.Feature) -> vcy.Feature:
    """Jump to the next exon following transcription direction instead of chromosome coordinate

    Arguments
    ---------
    feature: vcy.Feature
        An exonic feature

    Returns
    -------
    vcy.Feature:
        The next 3' exon following transcription direction

    Note
    ----
    It raises IndexError if the Feature is the 3'-most feature in the transcript model
    """
    if feature.transcript_model.chromstrand[-1] == "+":
        ix = feature.exin_no * 2
    else:
        ix = len(feature.transcript_model.list_features) - 2 * (feature.exin_no - 1) - 3
        if ix < 0:
            raise IndexError
    return feature.transcript_model.list_features[ix]


def closest_3prime(segment_match: vcy.SegmentMatch) -> int:
    """Calculate the closest distance walking on the transcript model to the 3'UTR

    Argument
    --------
    segment_match: vcy.SegmentMatch
        the sement from which 5' extremity calculate the distance

    Returns
    -------
    int:
        Distance in base pairs

    Note
    ----
    It skips all the introns but the one where the segment is mapping (if mapping is intronic)

    """
    dist23prime: int = 0
    if segment_match.feature.transcript_model.chromstrand[-1] == "+":
        if segment_match.maps_to_exon:
            curr_exon = segment_match.feature
            to_end_of_exon = curr_exon.end - segment_match.segment[0] + 1
        else:  # maps_to_intron
            curr_intron = segment_match.feature
            to_end_of_exon = curr_intron.end - segment_match.segment[0] + 1
            curr_exon = curr_intron.get_downstream_exon()
            to_end_of_exon += len(curr_exon)

        dist23prime += to_end_of_exon
        while True:
            try:
                curr_exon = jump_next_3p_exon(curr_exon)
                dist23prime += len(curr_exon)
            except IndexError:
                break
    else:  # "-" strand
        if segment_match.maps_to_exon:
            curr_exon = segment_match.feature
            to_end_of_exon = segment_match.segment[-1] - curr_exon.start + 1
        else:  # maps_to_intron
            curr_intron = segment_match.feature
            to_end_of_exon = segment_match.segment[-1] - curr_intron.start + 1
            curr_exon = curr_intron.get_upstream_exon()
            to_end_of_exon += len(curr_exon)

        dist23prime += to_end_of_exon
        while True:
            try:
                curr_exon = jump_next_3p_exon(curr_exon)
                dist23prime += len(curr_exon)
            except IndexError:
                break
    return dist23prime


def spliced_iter(segments_list: List[vcy.SegmentMatch], read_len: int=99) -> Iterable[vcy.SegmentMatch]:
    """Iterates over a list of segment matches grouping spliced segment in new special segment matches
    The goal i s to make the yielded output to be compatible with closest_3prime and further 3' mapping dist pipeline

    Arguments
    ---------
    segments_list: List[vcy.SegmentMatch]
        A list of segment matches objects

    read_len: int, default=99
        The lenght of a illumin read in the given technology

    Returns
    -------
    Iterable[vcy.SegmentMatch]:
        Either the original segment match or a fake Segment matche that will take the place of a set segment that are spliced

    Note
    ----
    This does not take in consideration all the corner cases it is very difficoult without keeping track of the splicing event
    """
    segments_list = list(segments_list)  # avoid popping way results from the original
    i = 0
    while len(segments_list):
        sm = segments_list.pop(0)
        if sm.is_spliced:
            sm_list = [sm]
            while segments_list[0].is_spliced:
                sm_list.append(segments_list.pop(0))
                if sum([s.segment[1] - s.segment[0] + 1 for s in sm_list]) + segments_list[0] > read_len:
                    break
            if len(segments_list) != 2:
                # just for safety let's ignore those counts otherwise we could make a mess
                continue
                
            if sm_list[0].feature.transcript_model.chromstrand[-1] == "+":
                if sm_list[-1].feature.kind == ord("i"):
                    i += len(sm_list)
                    yield vcy.SegmentMatch(segment=sm_list[0].segment,
                                           feature=sm_list[-1].feature)
                else:  # ord("e")
                    # in this way the distance should be calculated correctly
                    i += len(sm_list)
                    yield vcy.SegmentMatch(segment=(sm_list[-1].feature.start - (sm_list[0].segment[-1] - sm_list[0].segment[0]), -1),
                                           feature=sm_list[-1].feature)
            else:  # "-" strand
                if sm_list[0].feature.kind == ord("i"):
                    i += len(sm_list)
                    yield vcy.SegmentMatch(segment=sm_list[-1].segment,
                                           feature=sm_list[0].feature)
                else:  # ord("e")
                    i += len(sm_list)
                    yield vcy.SegmentMatch(segment=(-1, sm_list[-1].feature.end + (sm_list[0].segment[-1] - sm_list[0].segment[0]), -1),
                                           feature=sm_list[0].feature)
        else:
            i += 1
            yield sm
