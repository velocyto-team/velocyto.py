from typing import *
import velocyto as vcy
import numpy as np
import abc


class Logic(metaclass=abc.ABCMeta):
    """Base class from wich all the logics should inherit
    """

    def __init__(self) -> None:
        self.name = "Logic"

    @property
    def layers(self) -> List[str]:  # This should be overridden if a different set of layers is desired
        return ["spliced", "unspliced", "ambiguous"]

    @abc.abstractmethod  # This needs to be overridden
    def count(self, molitem: vcy.Molitem, cell_bcidx: int, spliced: np.ndarray,
              unspliced: np.ndarray, ambiguous: np.ndarray, geneid2ix: Dict[str, int]) -> None:
        """This methods will have to countain the core operatios of the logic to attribute a molecule to one of the cathergories
        
        Arguments
        ---------
        molitem: vcy.Molitem
            The :py:class:`vcy.Molitem` object to be considered by the logic
        cell_bcidx: int
            The cell index in the memory buffers below
        spliced: np.ndarray,
            The memory buffer that will be saved in the loom file after counting
        unspliced: np.ndarray
            The memory buffer that will be saved in the loom file after counting
        ambiguous: np.ndarray
            The memory buffer that will be saved in the loom file after counting
        geneid2ix: Dict[str, int]
            Dictionary containing the Acession of the genes mapping to its column index position
        
        
        Returns
        -------
        Nothing but it adds the molecule to the appropriate layer (or does not count at all)

        Note
        ----
        I need to generalize this to any set of layers: instead of spliced, unspliced  and ambiguous
        np.ndarray a dictionary {"name_layer": np.ndarray} should be passed
        """
        # NOTE I need to generalize this to any set of layers
        return


class Permissive10X(Logic):
    """Permissive logic for 10X Genomics chemistry

    This logic differs from the other 10x Logics because:
    - singletons if the fall in not validated introns are COUNTED UNSPLICED
    - singletons if the fall in validated introns are COUNTED UNSPLICED
    - non-singletons if are supported by not validated introns are COUNTED UNSPLICED
    - non-singletons if are supported by validated introns are COUNTED UNSPLICED
    """

    def __init__(self) -> None:
        self.name = "Permissive10X"

    @property
    def layers(self) -> List[str]:  # This should be overridden if a different set of layers is desired
        return ["spliced", "unspliced", "ambiguous"]

    def count(self, molitem: vcy.Molitem, cell_bcidx: int, spliced: np.ndarray,
              unspliced: np.ndarray, ambiguous: np.ndarray, geneid2ix: Dict[str, int]) -> None:

        # The hits are not compatible with any annotated transcript model
        if len(molitem.mappings_record) == 0:
            return
        # Compatible with one or more transcript models:
        else:
            # Check that there are not different possible genes ??
            if len(set(i.geneid for i in molitem.mappings_record.keys())) == 1:
                gene_check: Set[str] = set()
                
                has_onlyintron_model = 0
                has_only_span_exin_model = 1
                has_onlyintron_and_valid_model = 0
                has_valid_mixed_model = 0
                has_invalid_mixed_model = 0
                has_onlyexo_model = 0
                has_mixed_model = 0
                multi_gene = 0
                for transcript_model, segments_list in molitem.mappings_record.items():
                    gene_check.add(transcript_model.geneid)
                    if len(gene_check) > 1:
                        multi_gene = 1
                    has_introns = 0
                    has_exons = 0
                    has_exseg_with_spliced_flag = 0
                    has_validated_intron = 0
                    has_exin_intron_span = 0
                    has_non3prime = 0
                    for segment_match in segments_list:
                        if segment_match.maps_to_intron:
                            has_introns = 1
                            if segment_match.feature.is_validated:
                                has_validated_intron = 1
                                if segment_match.feature.end_overlaps_with_part_of(segment_match.segment):
                                    downstream_exon = segment_match.feature.get_downstream_exon()
                                    if downstream_exon.start_overlaps_with_part_of(segment_match.segment):
                                        has_exin_intron_span = 1
                                if segment_match.feature.start_overlaps_with_part_of(segment_match.segment):
                                    upstream_exon = segment_match.feature.get_upstream_exon()
                                    if upstream_exon.end_overlaps_with_part_of(segment_match.segment):
                                        has_exin_intron_span = 1
                        elif segment_match.maps_to_exon:
                            has_exons = 1
                            if not segment_match.feature.is_last_3prime:
                                has_non3prime = 1
                            if segment_match.is_spliced:
                                has_exseg_with_spliced_flag = 1
                    if has_validated_intron and not has_exons:
                        has_onlyintron_and_valid_model = 1
                    if has_introns and not has_exons:
                        has_onlyintron_model = 1
                    if has_exons and not has_introns:
                        has_onlyexo_model = 1
                    if has_exons and has_introns and not has_validated_intron and not has_exin_intron_span:
                        has_invalid_mixed_model = 1
                        has_mixed_model = 1
                    if has_exons and has_introns and has_validated_intron and not has_exin_intron_span:
                        has_valid_mixed_model = 1
                        has_mixed_model = 1
                    if not has_exin_intron_span:
                        has_only_span_exin_model = 0
                        
                if multi_gene:
                    # Many genes are compatible with the observation, do not count
                    return
                else:
                    if not len(molitem.mappings_record):
                        # No gene is compatible with the observation, do not count
                        return
                    else:
                        if has_onlyexo_model and not has_onlyintron_model and not has_mixed_model:
                            # More common situation, normal exonic read, count as spliced
                            gene_ix = geneid2ix[transcript_model.geneid]
                            spliced[gene_ix, cell_bcidx] += 1
                            return
                        if has_only_span_exin_model:
                            # All the compatible transcript models have spanning exon-intron boundaries, count unspliced
                            gene_ix = geneid2ix[transcript_model.geneid]
                            unspliced[gene_ix, cell_bcidx] += 1
                            return
                        if has_onlyintron_and_valid_model and not has_mixed_model and not has_onlyexo_model:
                            if len(segments_list) == 1:
                                # Singleton in validated intron
                                gene_ix = geneid2ix[transcript_model.geneid]
                                unspliced[gene_ix, cell_bcidx] += 1
                                return
                            else:
                                # Non-singleton in validated intron
                                gene_ix = geneid2ix[transcript_model.geneid]
                                unspliced[gene_ix, cell_bcidx] += 1
                                return
                        if has_onlyintron_model and not has_onlyintron_and_valid_model and not has_mixed_model and not has_onlyexo_model:
                            if len(segments_list) == 1:
                                # Singleton in non-validated intron
                                gene_ix = geneid2ix[transcript_model.geneid]
                                unspliced[gene_ix, cell_bcidx] += 1
                                return
                            else:
                                # Non-singleton in non-validated intron
                                gene_ix = geneid2ix[transcript_model.geneid]
                                unspliced[gene_ix, cell_bcidx] += 1
                                return
                        if has_invalid_mixed_model and not has_valid_mixed_model and not has_onlyintron_model and not has_onlyexo_model and not has_only_span_exin_model:
                            # Not validated and mapping to exon and introns, happens rarely in 10X / irrelevant. Count anyways
                            gene_ix = geneid2ix[transcript_model.geneid]
                            unspliced[gene_ix, cell_bcidx] += 1
                            return
                        if has_valid_mixed_model and not has_onlyintron_model and not has_onlyexo_model and not has_only_span_exin_model:
                            # Validated and mapping to exon and introns, happens rarely in 10X. Count as unspliced.
                            gene_ix = geneid2ix[transcript_model.geneid]
                            unspliced[gene_ix, cell_bcidx] += 1
                            return
                        if has_onlyintron_model and has_onlyexo_model and not has_mixed_model:
                            # Ambiguity among the transcript models compatible with the mapping, most common case! Count ambiguous
                            gene_ix = geneid2ix[transcript_model.geneid]
                            ambiguous[gene_ix, cell_bcidx] += 1
                            return
                        if has_onlyintron_model and not has_onlyexo_model and has_mixed_model:
                            # Ambiguity among the transcript models compatible with the mapping. Very rare. Count ambiguous
                            gene_ix = geneid2ix[transcript_model.geneid]
                            ambiguous[gene_ix, cell_bcidx] += 1
                            return
                        if not has_onlyintron_model and has_onlyexo_model and has_mixed_model:
                            # Ambiguity among the transcript models compatible with the mapping. Count ambiguous
                            gene_ix = geneid2ix[transcript_model.geneid]
                            ambiguous[gene_ix, cell_bcidx] += 1
                            return
                        if has_onlyintron_model and has_onlyexo_model and has_mixed_model:
                            # Ambiguity among the transcript models compatible with the mapping. Very rare. Count ambiguous
                            gene_ix = geneid2ix[transcript_model.geneid]
                            ambiguous[gene_ix, cell_bcidx] += 1
                            return


class ValidatedIntrons10X(Logic):
    """ValidatedIntrons logic for 10X Genomics chemistry
    
    This differs from the other 10x Logics because:
    - singletons if the fall in not validated introns are DISCARDED
    - singletons if the fall in validated introns are COUNTED UNSPLICED
    - non-singletons if are supported by not validated introns are DISCARDED
    - non-singletons if are supported by validated introns are COUNTED UNSPLICED

    """

    def __init__(self) -> None:
        self.name = "ValidatedIntrons10X"

    @property
    def layers(self) -> List[str]:  # This should be overridden if a different set of layers is desired
        return ["spliced", "unspliced", "ambiguous"]

    def count(self, molitem: vcy.Molitem, cell_bcidx: int, spliced: np.ndarray,
              unspliced: np.ndarray, ambiguous: np.ndarray, geneid2ix: Dict[str, int]) -> None:

        # The hits are not compatible with any annotated transcript model
        if len(molitem.mappings_record) == 0:
            return
        # Compatible with one or more transcript models:
        else:
            # Check that there are not different possible genes ??
            if len(set(i.geneid for i in molitem.mappings_record.keys())) == 1:
                gene_check: Set[str] = set()
                
                has_onlyintron_model = 0
                has_only_span_exin_model = 1
                has_onlyintron_and_valid_model = 0
                has_valid_mixed_model = 0
                has_invalid_mixed_model = 0
                has_onlyexo_model = 0
                has_mixed_model = 0
                multi_gene = 0
                for transcript_model, segments_list in molitem.mappings_record.items():
                    gene_check.add(transcript_model.geneid)
                    if len(gene_check) > 1:
                        multi_gene = 1
                    has_introns = 0
                    has_exons = 0
                    has_exseg_with_spliced_flag = 0
                    has_validated_intron = 0
                    has_exin_intron_span = 0
                    has_non3prime = 0
                    for segment_match in segments_list:
                        if segment_match.maps_to_intron:
                            has_introns = 1
                            if segment_match.feature.is_validated:
                                has_validated_intron = 1
                                if segment_match.feature.end_overlaps_with_part_of(segment_match.segment):
                                    downstream_exon = segment_match.feature.get_downstream_exon()
                                    if downstream_exon.start_overlaps_with_part_of(segment_match.segment):
                                        has_exin_intron_span = 1
                                if segment_match.feature.start_overlaps_with_part_of(segment_match.segment):
                                    upstream_exon = segment_match.feature.get_upstream_exon()
                                    if upstream_exon.end_overlaps_with_part_of(segment_match.segment):
                                        has_exin_intron_span = 1
                        elif segment_match.maps_to_exon:
                            has_exons = 1
                            if not segment_match.feature.is_last_3prime:
                                has_non3prime = 1
                            if segment_match.is_spliced:
                                has_exseg_with_spliced_flag = 1
                    if has_validated_intron and not has_exons:
                        has_onlyintron_and_valid_model = 1
                    if has_introns and not has_exons:
                        has_onlyintron_model = 1
                    if has_exons and not has_introns:
                        has_onlyexo_model = 1
                    if has_exons and has_introns and not has_validated_intron and not has_exin_intron_span:
                        has_invalid_mixed_model = 1
                        has_mixed_model = 1
                    if has_exons and has_introns and has_validated_intron and not has_exin_intron_span:
                        has_valid_mixed_model = 1
                        has_mixed_model = 1
                    if not has_exin_intron_span:
                        has_only_span_exin_model = 0
                        
                if multi_gene:
                    # many genes are compatible with the observation, do not count
                    return
                else:
                    if not len(molitem.mappings_record):
                        # no gene is compatible with the observation, do not count
                        return
                    else:
                        if has_onlyexo_model and not has_onlyintron_model and not has_mixed_model:
                            # More common situation, normal exonic read, count as spliced
                            gene_ix = geneid2ix[transcript_model.geneid]
                            spliced[gene_ix, cell_bcidx] += 1
                            return
                        if has_only_span_exin_model:
                            # All the compatible transcript models have spanning exon-intron boundaries, count unspliced
                            gene_ix = geneid2ix[transcript_model.geneid]
                            unspliced[gene_ix, cell_bcidx] += 1
                            return
                        if has_onlyintron_and_valid_model and not has_mixed_model and not has_onlyexo_model:
                            if len(segments_list) == 1:
                                # Singleton in validated intron
                                gene_ix = geneid2ix[transcript_model.geneid]
                                unspliced[gene_ix, cell_bcidx] += 1
                                return
                            else:
                                # Non singleton in validated intron
                                gene_ix = geneid2ix[transcript_model.geneid]
                                unspliced[gene_ix, cell_bcidx] += 1
                                return
                        if has_onlyintron_model and not has_onlyintron_and_valid_model and not has_mixed_model and not has_onlyexo_model:
                            if len(segments_list) == 1:
                                # Singleton in non-validated intron
                                return
                            else:
                                # Non-singleton in non-validated intron
                                return
                        if has_invalid_mixed_model and not has_valid_mixed_model and not has_onlyintron_model and not has_onlyexo_model and not has_only_span_exin_model:
                            # Not validated and mapping to exon and introns, happens rarely in 10X / irrelevant.
                            return
                        if has_valid_mixed_model and not has_onlyintron_model and not has_onlyexo_model and not has_only_span_exin_model:
                            # Validated and mapping to exon and introns, happens rarely in 10X. Count as unspliced.
                            gene_ix = geneid2ix[transcript_model.geneid]
                            unspliced[gene_ix, cell_bcidx] += 1
                            return
                        if has_onlyintron_model and has_onlyexo_model and not has_mixed_model:
                            # Ambiguity among the transcript models compatible with the mapping, most common case! Count ambiguous
                            gene_ix = geneid2ix[transcript_model.geneid]
                            ambiguous[gene_ix, cell_bcidx] += 1
                            return
                        if has_onlyintron_model and not has_onlyexo_model and has_mixed_model:
                            # Ambiguity among the transcript models compatible with the mapping. Very rare.
                            gene_ix = geneid2ix[transcript_model.geneid]
                            ambiguous[gene_ix, cell_bcidx] += 1
                            return
                        if not has_onlyintron_model and has_onlyexo_model and has_mixed_model:
                            # Ambiguity among the transcript models compatible with the mapping. Count ambiguous
                            gene_ix = geneid2ix[transcript_model.geneid]
                            ambiguous[gene_ix, cell_bcidx] += 1
                            return
                        if has_onlyintron_model and has_onlyexo_model and has_mixed_model:
                            # Ambiguity among the transcript models compatible with the mapping. Very rare.
                            gene_ix = geneid2ix[transcript_model.geneid]
                            ambiguous[gene_ix, cell_bcidx] += 1
                            return


class Stricter10X(Logic):
    """Stricter logic for 10X Genomics chemistry

    This differ from the other 10x Logics because:
    - singletons if the fall in not validated introns are DISCARDED
    - singletons if the fall in validated introns are DISCARDED
    - non-singletons if are supported by not validated introns are DISCARDED
    - non-singletons if are supported by validated introns are COUNTED UNSPLICED

    """

    def __init__(self) -> None:
        self.name = "Stricter10X"

    @property
    def layers(self) -> List[str]:  # This should be overridden if a different set of layers is desired
        return ["spliced", "unspliced", "ambiguous"]

    def count(self, molitem: vcy.Molitem, cell_bcidx: int, spliced: np.ndarray,
              unspliced: np.ndarray, ambiguous: np.ndarray, geneid2ix: Dict[str, int]) -> None:

        # The hits are not compatible with any annotated transcript model
        if len(molitem.mappings_record) == 0:
            return
        # Compatible with one or more transcript models:
        else:
            # Check that there are not different possible genes ??
            if len(set(i.geneid for i in molitem.mappings_record.keys())) == 1:
                gene_check: Set[str] = set()
                
                has_onlyintron_model = 0
                has_only_span_exin_model = 1
                has_onlyintron_and_valid_model = 0
                has_valid_mixed_model = 0
                has_invalid_mixed_model = 0
                has_onlyexo_model = 0
                has_mixed_model = 0
                multi_gene = 0
                for transcript_model, segments_list in molitem.mappings_record.items():
                    gene_check.add(transcript_model.geneid)
                    if len(gene_check) > 1:
                        multi_gene = 1
                    has_introns = 0
                    has_exons = 0
                    has_exseg_with_spliced_flag = 0
                    has_validated_intron = 0
                    has_exin_intron_span = 0
                    has_non3prime = 0
                    for segment_match in segments_list:
                        if segment_match.maps_to_intron:
                            has_introns = 1
                            if segment_match.feature.is_validated:
                                has_validated_intron = 1
                                if segment_match.feature.end_overlaps_with_part_of(segment_match.segment):
                                    downstream_exon = segment_match.feature.get_downstream_exon()
                                    if downstream_exon.start_overlaps_with_part_of(segment_match.segment):
                                        has_exin_intron_span = 1
                                if segment_match.feature.start_overlaps_with_part_of(segment_match.segment):
                                    upstream_exon = segment_match.feature.get_upstream_exon()
                                    if upstream_exon.end_overlaps_with_part_of(segment_match.segment):
                                        has_exin_intron_span = 1
                        elif segment_match.maps_to_exon:
                            has_exons = 1
                            if not segment_match.feature.is_last_3prime:
                                has_non3prime = 1
                            if segment_match.is_spliced:
                                has_exseg_with_spliced_flag = 1
                    if has_validated_intron and not has_exons:
                        has_onlyintron_and_valid_model = 1
                    if has_introns and not has_exons:
                        has_onlyintron_model = 1
                    if has_exons and not has_introns:
                        has_onlyexo_model = 1
                    if has_exons and has_introns and not has_validated_intron and not has_exin_intron_span:
                        has_invalid_mixed_model = 1
                        has_mixed_model = 1
                    if has_exons and has_introns and has_validated_intron and not has_exin_intron_span:
                        has_valid_mixed_model = 1
                        has_mixed_model = 1
                    if not has_exin_intron_span:
                        has_only_span_exin_model = 0
                        
                if multi_gene:
                    # many genes are compatible with the observation, do not count
                    return
                else:
                    if not len(molitem.mappings_record):
                        # no gene is compatible with the observation, do not count
                        return
                    else:
                        if has_onlyexo_model and not has_onlyintron_model and not has_mixed_model:
                            # More common situation, normal exonic read, count as spliced
                            gene_ix = geneid2ix[transcript_model.geneid]
                            spliced[gene_ix, cell_bcidx] += 1
                            return
                        if has_only_span_exin_model:
                            # All the compatible transcript models have spanning exon-intron boundaries, count unspliced
                            gene_ix = geneid2ix[transcript_model.geneid]
                            unspliced[gene_ix, cell_bcidx] += 1
                            return
                        if has_onlyintron_and_valid_model and not has_mixed_model and not has_onlyexo_model:
                            if len(segments_list) == 1:
                                # Singleton in validated intron, do not count
                                return
                            else:
                                # Non singleton in validated intron
                                gene_ix = geneid2ix[transcript_model.geneid]
                                unspliced[gene_ix, cell_bcidx] += 1
                                return
                        if has_onlyintron_model and not has_onlyintron_and_valid_model and not has_mixed_model and not has_onlyexo_model:
                            if len(segments_list) == 1:
                                # Singleton in non-validated intron
                                return
                            else:
                                # Non-singleton in non-validated intron
                                return
                        if has_invalid_mixed_model and not has_valid_mixed_model and not has_onlyintron_model and not has_onlyexo_model and not has_only_span_exin_model:
                            # Not validated and mapping to exon and introns, happens rarely in 10X / irrelevant.
                            return
                        if has_valid_mixed_model and not has_onlyintron_model and not has_onlyexo_model and not has_only_span_exin_model:
                            # Validated and mapping to exon and introns, happens rarely in 10X. Count as unspliced.
                            gene_ix = geneid2ix[transcript_model.geneid]
                            unspliced[gene_ix, cell_bcidx] += 1
                            return
                        if has_onlyintron_model and has_onlyexo_model and not has_mixed_model:
                            # Ambiguity among the transcript models compatible with the mapping, most common case! Count ambiguous
                            gene_ix = geneid2ix[transcript_model.geneid]
                            ambiguous[gene_ix, cell_bcidx] += 1
                            return
                        if has_onlyintron_model and not has_onlyexo_model and has_mixed_model:
                            # Ambiguity among the transcript models compatible with the mapping. Very rare.
                            gene_ix = geneid2ix[transcript_model.geneid]
                            ambiguous[gene_ix, cell_bcidx] += 1
                            return
                        if not has_onlyintron_model and has_onlyexo_model and has_mixed_model:
                            # Ambiguity among the transcript models compatible with the mapping. Count ambiguous
                            gene_ix = geneid2ix[transcript_model.geneid]
                            ambiguous[gene_ix, cell_bcidx] += 1
                            return
                        if has_onlyintron_model and has_onlyexo_model and has_mixed_model:
                            # Ambiguity among the transcript models compatible with the mapping. Very rare.
                            gene_ix = geneid2ix[transcript_model.geneid]
                            ambiguous[gene_ix, cell_bcidx] += 1
                            return


class ObservedSpanning10X(Logic):
    """ObservedSpanning10X logic for 10X Genomics chemistry

    This differ from the other 10x Logics because:
    - singletons if the fall in not validated introns are DISCARDED
    - singletons if the fall in validated introns are DISCARDED
    - non-singletons if are supported by not validated introns are DISCARDED
    - non-singletons if are supported by validated introns are DISCARDED
    - Therefore only the observed intron spanning reads are counted as UNSPLICED
    """

    def __init__(self) -> None:
        self.name = "ObservedSpanning10X"

    @property
    def layers(self) -> List[str]:  # This should be overridden if a different set of layers is desired
        return ["spliced", "unspliced", "ambiguous"]

    def count(self, molitem: vcy.Molitem, cell_bcidx: int, spliced: np.ndarray,
              unspliced: np.ndarray, ambiguous: np.ndarray, geneid2ix: Dict[str, int]) -> None:

        # The hits are not compatible with any annotated transcript model
        if len(molitem.mappings_record) == 0:
            return
        # Compatible with one or more transcript models:
        else:
            # Check that there are not different possible genes ??
            if len(set(i.geneid for i in molitem.mappings_record.keys())) == 1:
                gene_check: Set[str] = set()
                
                has_onlyintron_model = 0
                has_only_span_exin_model = 1
                has_onlyintron_and_valid_model = 0
                has_valid_mixed_model = 0
                has_invalid_mixed_model = 0
                has_onlyexo_model = 0
                has_mixed_model = 0
                multi_gene = 0
                for transcript_model, segments_list in molitem.mappings_record.items():
                    gene_check.add(transcript_model.geneid)
                    if len(gene_check) > 1:
                        multi_gene = 1
                    has_introns = 0
                    has_exons = 0
                    has_exseg_with_spliced_flag = 0
                    has_validated_intron = 0
                    has_exin_intron_span = 0
                    has_non3prime = 0
                    for segment_match in segments_list:
                        if segment_match.maps_to_intron:
                            has_introns = 1
                            if segment_match.feature.is_validated:
                                has_validated_intron = 1
                                if segment_match.feature.end_overlaps_with_part_of(segment_match.segment):
                                    downstream_exon = segment_match.feature.get_downstream_exon()
                                    if downstream_exon.start_overlaps_with_part_of(segment_match.segment):
                                        has_exin_intron_span = 1
                                if segment_match.feature.start_overlaps_with_part_of(segment_match.segment):
                                    upstream_exon = segment_match.feature.get_upstream_exon()
                                    if upstream_exon.end_overlaps_with_part_of(segment_match.segment):
                                        has_exin_intron_span = 1
                        elif segment_match.maps_to_exon:
                            has_exons = 1
                            if not segment_match.feature.is_last_3prime:
                                has_non3prime = 1
                            if segment_match.is_spliced:
                                has_exseg_with_spliced_flag = 1
                    if has_validated_intron and not has_exons:
                        has_onlyintron_and_valid_model = 1
                    if has_introns and not has_exons:
                        has_onlyintron_model = 1
                    if has_exons and not has_introns:
                        has_onlyexo_model = 1
                    if has_exons and has_introns and not has_validated_intron and not has_exin_intron_span:
                        has_invalid_mixed_model = 1
                        has_mixed_model = 1
                    if has_exons and has_introns and has_validated_intron and not has_exin_intron_span:
                        has_valid_mixed_model = 1
                        has_mixed_model = 1
                    if not has_exin_intron_span:
                        has_only_span_exin_model = 0
                        
                if multi_gene:
                    # Many genes are compatible with the observation, do not count
                    return
                else:
                    if not len(molitem.mappings_record):
                        # No gene is compatible with the observation, do not count
                        return
                    else:
                        if has_onlyexo_model and not has_onlyintron_model and not has_mixed_model:
                            # More common situation, normal exonic read, count as spliced
                            gene_ix = geneid2ix[transcript_model.geneid]
                            spliced[gene_ix, cell_bcidx] += 1
                            return
                        if has_only_span_exin_model:
                            # All the compatible transcript models have spanning exon-intron boundaries, count unspliced
                            gene_ix = geneid2ix[transcript_model.geneid]
                            unspliced[gene_ix, cell_bcidx] += 1
                            return
                        if has_onlyintron_and_valid_model and not has_mixed_model and not has_onlyexo_model:
                            if len(segments_list) == 1:
                                # Singleton in validated intron, do not count
                                return
                            else:
                                # Non-singleton in validated intron, do not count
                                return
                        if has_onlyintron_model and not has_onlyintron_and_valid_model and not has_mixed_model and not has_onlyexo_model:
                            if len(segments_list) == 1:
                                # Singleton in non-validated intron
                                return
                            else:
                                # Non-singleton in non-validated intron
                                return
                        if has_invalid_mixed_model and not has_valid_mixed_model and not has_onlyintron_model and not has_onlyexo_model and not has_only_span_exin_model:
                            # Not validated and mapping to exon and introns, happens rarely in 10X / irrelevant.
                            return
                        if has_valid_mixed_model and not has_onlyintron_model and not has_onlyexo_model and not has_only_span_exin_model:
                            # Validated and mapping to exon and introns, happens rarely in 10X. Count as unspliced.
                            gene_ix = geneid2ix[transcript_model.geneid]
                            unspliced[gene_ix, cell_bcidx] += 1
                            return
                        if has_onlyintron_model and has_onlyexo_model and not has_mixed_model:
                            # Ambiguity among the transcript models compatible with the mapping, most common case! Count ambiguous
                            gene_ix = geneid2ix[transcript_model.geneid]
                            ambiguous[gene_ix, cell_bcidx] += 1
                            return
                        if has_onlyintron_model and not has_onlyexo_model and has_mixed_model:
                            # Ambiguity among the transcript models compatible with the mapping. Very rare.
                            gene_ix = geneid2ix[transcript_model.geneid]
                            ambiguous[gene_ix, cell_bcidx] += 1
                            return
                        if not has_onlyintron_model and has_onlyexo_model and has_mixed_model:
                            # Ambiguity among the transcript models compatible with the mapping. Count ambiguous
                            gene_ix = geneid2ix[transcript_model.geneid]
                            ambiguous[gene_ix, cell_bcidx] += 1
                            return
                        if has_onlyintron_model and has_onlyexo_model and has_mixed_model:
                            # Ambiguity among the transcript models compatible with the mapping. Very rare.
                            gene_ix = geneid2ix[transcript_model.geneid]
                            ambiguous[gene_ix, cell_bcidx] += 1
                            return


Default: type = ValidatedIntrons10X
