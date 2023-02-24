import abc
from typing import Union

import numpy as np

from .molitem import Molitem


class Logic(metaclass=abc.ABCMeta):
    """Base class from wich all the logics should inherit"""

    def __init__(self) -> None:
        self.name = "Logic"

    @property
    def layers(
        self,
    ) -> list[str]:  # This should be overridden if a different set of layers is desired
        return []

    @property
    def stranded(self) -> bool:
        return True

    @property
    def perform_validation_markup(self) -> bool:
        return True

    @property
    def accept_discordant(self) -> bool:
        return False

    @abc.abstractmethod  # This needs to be overridden
    def count(
        self,
        molitem: Molitem,
        cell_bcidx: int,
        dict_layers_columns: dict[str, np.ndarray],
        geneid2ix: dict[str, int],
    ) -> Union[None, int]:
        """This methods will have to countain the core operations of the logic to attribute a molecule to one of the cathergories

        Arguments
        ---------
        molitem: Molitem
            The :py:class:`Molitem` object to be considered by the logic
        cell_bcidx: int
            The cell index in the memory buffers below
        dict_layers_columns: dict[str, np.ndarray]
            A dictionary mapping the name of a layer with the memory buffer that will be saved in the loom file after counting
        geneid2ix: dict[str, int]
            dictionary containing the Acession of the genes mapping to its column index position

        Returns
        -------
        Nothing but it adds the molecule to the appropriate layer (or does not count at all)

        """
        # NOTE I need to generalize this to any set of layers
        return None


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
    def layers(
        self,
    ) -> list[str]:  # This should be overridden if a different set of layers is desired
        return ["spliced", "unspliced", "ambiguous"]

    @property
    def stranded(self) -> bool:
        return True

    @property
    def perform_validation_markup(self) -> bool:
        return True

    @property
    def accept_discordant(self) -> bool:
        return False

    def count(
        self,
        molitem: Molitem,
        cell_bcidx: int,
        dict_layers_columns: dict[str, np.ndarray],
        geneid2ix: dict[str, int],
    ) -> int:
        if len(molitem.mappings_record) == 0:
            return 2
        # NOTE This can be simplified qyuite a bit, without loss of acuracy!
        # The hits are not compatible with any annotated transcript model
        spliced = dict_layers_columns["spliced"]  # .toarray()
        unspliced = dict_layers_columns["unspliced"]  # .toarray()
        ambiguous = dict_layers_columns["ambiguous"]  # .toarray()
        if len({i.geneid for i in molitem.mappings_record.keys()}) != 1:
            return 3
        gene_check: set[str] = set()

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
            # has_exseg_with_spliced_flag = 0
            has_validated_intron = 0
            has_exin_intron_span = 0
            # has_non3prime = 0
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
                    # if not segment_match.feature.is_last_3prime:
                    #     has_non3prime = 1
                    # if segment_match.is_spliced:
                    #     has_exseg_with_spliced_flag = 1
            if has_validated_intron and not has_exons:
                has_onlyintron_and_valid_model = 1
            if has_introns and not has_exons:
                has_onlyintron_model = 1
            if has_exons:
                if not has_introns:
                    has_onlyexo_model = 1
                if has_introns:
                    if not has_validated_intron and not has_exin_intron_span:
                        has_invalid_mixed_model = 1
                        has_mixed_model = 1
                    if has_validated_intron and not has_exin_intron_span:
                        has_valid_mixed_model = 1
                        has_mixed_model = 1
            if not has_exin_intron_span:
                has_only_span_exin_model = 0

        if multi_gene:
            # Many genes are compatible with the observation, do not count
            return 1
        if not len(molitem.mappings_record):
            # No gene is compatible with the observation, do not count
            return 2
        if has_onlyexo_model and not has_onlyintron_model and not has_mixed_model:
            return self._extracted_from_count_(
                geneid2ix, transcript_model, spliced, cell_bcidx
            )
        if has_only_span_exin_model:
            return self._extracted_from_count_(
                geneid2ix, transcript_model, unspliced, cell_bcidx
            )
        if has_onlyintron_and_valid_model and not has_mixed_model and not has_onlyexo_model:
            return self._extracted_from_count_(
                geneid2ix, transcript_model, unspliced, cell_bcidx
            )
        if (
                            has_onlyintron_model
                            and not has_onlyintron_and_valid_model
                            and not has_mixed_model
                            and not has_onlyexo_model
                        ):
            return self._extracted_from_count_(
                geneid2ix, transcript_model, unspliced, cell_bcidx
            )
        if (
            has_invalid_mixed_model
            and not has_valid_mixed_model
            and not has_onlyintron_model
            and not has_onlyexo_model
        ):
            return self._extracted_from_count_(
                geneid2ix, transcript_model, unspliced, cell_bcidx
            )
        if (
            has_valid_mixed_model
            and not has_onlyintron_model
            and not has_onlyexo_model
        ):
            return self._extracted_from_count_(
                geneid2ix, transcript_model, unspliced, cell_bcidx
            )
        if has_onlyintron_model:
            if has_onlyexo_model and not has_mixed_model:
                return self._extracted_from_count_(
                    geneid2ix, transcript_model, ambiguous, cell_bcidx
                )
            if not has_onlyexo_model:
                return self._extracted_from_count_(
                    geneid2ix, transcript_model, unspliced, cell_bcidx
                )
        if not has_onlyintron_model and has_onlyexo_model:
            return self._extracted_from_count_(
                geneid2ix, transcript_model, ambiguous, cell_bcidx
            )
        if has_onlyintron_model:
            return self._extracted_from_count_(
                geneid2ix, transcript_model, ambiguous, cell_bcidx
            )
        return 4

    # TODO Rename this here and in `count`
    def _extracted_from_count_(self, geneid2ix, transcript_model, arg2, cell_bcidx):
        # More common situation, normal exonic read, count as spliced
        gene_ix = geneid2ix[transcript_model.geneid]
        arg2[gene_ix, cell_bcidx] += 1
        return 0


class Intermediate10X(Logic):
    """ValidatedIntrons logic for 10X Genomics chemistry

    This differs from the other 10x Logics because:
    - singletons if the fall in not validated introns are DISCARDED
    - singletons if the fall in validated introns are COUNTED UNSPLICED
    - non-singletons if are supported by not validated introns are COUNTED UNSPLICED
    - non-singletons if are supported by validated introns are COUNTED UNSPLICED

    """

    def __init__(self) -> None:
        self.name = "Intermediate10X"

    @property
    def layers(
        self,
    ) -> list[str]:  # This should be overridden if a different set of layers is desired
        return ["spliced", "unspliced", "ambiguous"]

    @property
    def stranded(self) -> bool:
        return True

    @property
    def perform_validation_markup(self) -> bool:
        return True

    @property
    def accept_discordant(self) -> bool:
        return False

    def count(
        self,
        molitem: Molitem,
        cell_bcidx: int,
        dict_layers_columns: dict[str, np.ndarray],
        geneid2ix: dict[str, int],
    ) -> None:
        if len(molitem.mappings_record) == 0:
            return
        # NOTE This can be simplified qyuite a bit, without loss of acuracy!
        # The hits are not compatible with any annotated transcript model
        spliced = dict_layers_columns["spliced"]
        unspliced = dict_layers_columns["unspliced"]
        ambiguous = dict_layers_columns["ambiguous"]
            # Check that there are not different possible genes ??
        if len({i.geneid for i in molitem.mappings_record.keys()}) == 1:
            gene_check: set[str] = set()

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
                # has_exseg_with_spliced_flag = 0
                has_validated_intron = 0
                has_exin_intron_span = 0
                # has_non3prime = 0
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
                        # if not segment_match.feature.is_last_3prime:
                        #     has_non3prime = 1
                        # if segment_match.is_spliced:
                        #     has_exseg_with_spliced_flag = 1
                if has_validated_intron and not has_exons:
                    has_onlyintron_and_valid_model = 1
                if has_introns and not has_exons:
                    has_onlyintron_model = 1
                if has_exons:
                    if not has_introns:
                        has_onlyexo_model = 1
                    if has_introns:
                        if not has_validated_intron and not has_exin_intron_span:
                            has_invalid_mixed_model = 1
                            has_mixed_model = 1
                        if has_validated_intron and not has_exin_intron_span:
                            has_valid_mixed_model = 1
                            has_mixed_model = 1
                if not has_exin_intron_span:
                    has_only_span_exin_model = 0

            if not multi_gene and not len(molitem.mappings_record) or multi_gene:
                # no gene is compatible with the observation, do not count
                return
            if has_onlyexo_model and not has_onlyintron_model and not has_mixed_model:
                return self._extracted_from_count_(
                    geneid2ix, transcript_model, spliced, cell_bcidx
                )
            if has_only_span_exin_model:
                return self._extracted_from_count_(
                    geneid2ix, transcript_model, unspliced, cell_bcidx
                )
            if has_onlyintron_and_valid_model and not has_mixed_model and not has_onlyexo_model:
                return self._extracted_from_count_(
                    geneid2ix, transcript_model, unspliced, cell_bcidx
                )
            if (
                            has_onlyintron_model
                            and not has_onlyintron_and_valid_model
                            and not has_mixed_model
                            and not has_onlyexo_model
                        ):
                if len(segments_list) != 1:
                    # Non-singleton in non-validated intron
                    gene_ix = geneid2ix[transcript_model.geneid]
                    unspliced[gene_ix, cell_bcidx] += 1
                # Singleton in non-validated intron
                return
            if (
                has_invalid_mixed_model
                and not has_valid_mixed_model
                and not has_onlyintron_model
                and not has_onlyexo_model
            ):
                # Not validated and mapping to exon and introns, happens rarely in 10X / irrelevant.
                return
            if (
                has_valid_mixed_model
                and not has_onlyintron_model
                and not has_onlyexo_model
            ):
                return self._extracted_from_count_(
                    geneid2ix, transcript_model, unspliced, cell_bcidx
                )
            if has_onlyintron_model:
                if has_onlyexo_model and not has_mixed_model:
                    return self._extracted_from_count_(
                        geneid2ix, transcript_model, ambiguous, cell_bcidx
                    )
                if not has_onlyexo_model:
                    return self._extracted_from_count_(
                        geneid2ix, transcript_model, ambiguous, cell_bcidx
                    )
            if not has_onlyintron_model and has_onlyexo_model:
                return self._extracted_from_count_(
                    geneid2ix, transcript_model, ambiguous, cell_bcidx
                )
            if has_onlyintron_model:
                return self._extracted_from_count_(
                    geneid2ix, transcript_model, ambiguous, cell_bcidx
                )

    # TODO Rename this here and in `count`
    def _extracted_from_count_(self, geneid2ix, transcript_model, arg2, cell_bcidx):
        # More common situation, normal exonic read, count as spliced
        gene_ix = geneid2ix[transcript_model.geneid]
        arg2[gene_ix, cell_bcidx] += 1
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
    def layers(
        self,
    ) -> list[str]:  # This should be overridden if a different set of layers is desired
        return ["spliced", "unspliced", "ambiguous"]

    @property
    def stranded(self) -> bool:
        return True

    @property
    def perform_validation_markup(self) -> bool:
        return True

    @property
    def accept_discordant(self) -> bool:
        return False

    def count(
        self,
        molitem: Molitem,
        cell_bcidx: int,
        dict_layers_columns: dict[str, np.ndarray],
        geneid2ix: dict[str, int],
    ) -> None:
        if len(molitem.mappings_record) == 0:
            return
        # NOTE This can be simplified qyuite a bit, without loss of acuracy!
        # The hits are not compatible with any annotated transcript model
        spliced = dict_layers_columns["spliced"]
        unspliced = dict_layers_columns["unspliced"]
        ambiguous = dict_layers_columns["ambiguous"]
        if len({i.geneid for i in molitem.mappings_record.keys()}) == 1:
            gene_check: set[str] = set()

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
                # has_exseg_with_spliced_flag = 0
                has_validated_intron = 0
                has_exin_intron_span = 0
                # has_non3prime = 0
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
                        # if not segment_match.feature.is_last_3prime:
                        #     has_non3prime = 1
                        # if segment_match.is_spliced:
                        #     has_exseg_with_spliced_flag = 1
                if has_validated_intron and not has_exons:
                    has_onlyintron_and_valid_model = 1
                if has_introns and not has_exons:
                    has_onlyintron_model = 1
                if has_exons:
                    if not has_introns:
                        has_onlyexo_model = 1
                    if has_introns:
                        if not has_validated_intron and not has_exin_intron_span:
                            has_invalid_mixed_model = 1
                            has_mixed_model = 1
                        if has_validated_intron and not has_exin_intron_span:
                            has_valid_mixed_model = 1
                            has_mixed_model = 1
                if not has_exin_intron_span:
                    has_only_span_exin_model = 0

            if not multi_gene and not len(molitem.mappings_record) or multi_gene:
                # no gene is compatible with the observation, do not count
                return
            if has_onlyexo_model and not has_onlyintron_model and not has_mixed_model:
                return self._extracted_from_count_78(
                    geneid2ix, transcript_model, spliced, cell_bcidx
                )
            if has_only_span_exin_model:
                return self._extracted_from_count_78(
                    geneid2ix, transcript_model, unspliced, cell_bcidx
                )
            if has_onlyintron_and_valid_model and not has_mixed_model and not has_onlyexo_model:
                return self._extracted_from_count_78(
                    geneid2ix, transcript_model, unspliced, cell_bcidx
                )
            if (
                            has_onlyintron_model
                            and not has_onlyintron_and_valid_model
                            and not has_mixed_model
                            and not has_onlyexo_model
                        ):
                # Singleton in non-validated intron
                return
            if (
                has_invalid_mixed_model
                and not has_valid_mixed_model
                and not has_onlyintron_model
                and not has_onlyexo_model
            ):
                # Not validated and mapping to exon and introns, happens rarely in 10X / irrelevant.
                return
            if (
                has_valid_mixed_model
                and not has_onlyintron_model
                and not has_onlyexo_model
            ):
                return self._extracted_from_count_78(
                    geneid2ix, transcript_model, unspliced, cell_bcidx
                )
            if has_onlyintron_model:
                if has_onlyexo_model and not has_mixed_model:
                    return self._extracted_from_count_78(
                        geneid2ix, transcript_model, ambiguous, cell_bcidx
                    )
                if not has_onlyexo_model:
                    return self._extracted_from_count_78(
                        geneid2ix, transcript_model, ambiguous, cell_bcidx
                    )
            if not has_onlyintron_model and has_onlyexo_model:
                return self._extracted_from_count_78(
                    geneid2ix, transcript_model, ambiguous, cell_bcidx
                )
            if has_onlyintron_model:
                return self._extracted_from_count_78(
                    geneid2ix, transcript_model, ambiguous, cell_bcidx
                )

    # TODO Rename this here and in `count`
    def _extracted_from_count_78(self, geneid2ix, transcript_model, arg2, cell_bcidx):
        # More common situation, normal exonic read, count as spliced
        gene_ix = geneid2ix[transcript_model.geneid]
        arg2[gene_ix, cell_bcidx] += 1
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
    def layers(
        self,
    ) -> list[str]:  # This should be overridden if a different set of layers is desired
        return ["spliced", "unspliced", "ambiguous"]

    @property
    def stranded(self) -> bool:
        return True

    @property
    def perform_validation_markup(self) -> bool:
        return True

    def count(
        self,
        molitem: Molitem,
        cell_bcidx: int,
        dict_layers_columns: dict[str, np.ndarray],
        geneid2ix: dict[str, int],
    ) -> None:
        if len(molitem.mappings_record) == 0:
            return
        # NOTE This can be simplified qyuite a bit, without loss of acuracy!
        # The hits are not compatible with any annotated transcript model
        spliced = dict_layers_columns["spliced"]
        unspliced = dict_layers_columns["unspliced"]
        ambiguous = dict_layers_columns["ambiguous"]
            # Check that there are not different possible genes ??
        if len({i.geneid for i in molitem.mappings_record.keys()}) == 1:
            gene_check: set[str] = set()

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
                # has_exseg_with_spliced_flag = 0
                has_validated_intron = 0
                has_exin_intron_span = 0
                # has_non3prime = 0
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
                        # if not segment_match.feature.is_last_3prime:
                        #     has_non3prime = 1
                        # if segment_match.is_spliced:
                        #     has_exseg_with_spliced_flag = 1
                if has_validated_intron and not has_exons:
                    has_onlyintron_and_valid_model = 1
                if has_introns and not has_exons:
                    has_onlyintron_model = 1
                if has_exons:
                    if not has_introns:
                        has_onlyexo_model = 1
                    if has_introns:
                        if not has_validated_intron and not has_exin_intron_span:
                            has_invalid_mixed_model = 1
                            has_mixed_model = 1
                        if has_validated_intron and not has_exin_intron_span:
                            has_valid_mixed_model = 1
                            has_mixed_model = 1
                if not has_exin_intron_span:
                    has_only_span_exin_model = 0

            if not multi_gene and not len(molitem.mappings_record) or multi_gene:
                # no gene is compatible with the observation, do not count
                return
            if has_onlyexo_model and not has_onlyintron_model and not has_mixed_model:
                return self._extracted_from_count_(
                    geneid2ix, transcript_model, spliced, cell_bcidx
                )
            if has_only_span_exin_model:
                return self._extracted_from_count_(
                    geneid2ix, transcript_model, unspliced, cell_bcidx
                )
            if has_onlyintron_and_valid_model and not has_mixed_model and not has_onlyexo_model:
                if len(segments_list) != 1:
                    # Non singleton in validated intron
                    gene_ix = geneid2ix[transcript_model.geneid]
                    unspliced[gene_ix, cell_bcidx] += 1
                # Singleton in validated intron, do not count
                return
            if (
                            has_onlyintron_model
                            and not has_onlyintron_and_valid_model
                            and not has_mixed_model
                            and not has_onlyexo_model
                        ):
                # Singleton in non-validated intron
                return
            if (
                has_invalid_mixed_model
                and not has_valid_mixed_model
                and not has_onlyintron_model
                and not has_onlyexo_model
            ):
                # Not validated and mapping to exon and introns, happens rarely in 10X / irrelevant.
                return
            if (
                has_valid_mixed_model
                and not has_onlyintron_model
                and not has_onlyexo_model
            ):
                return self._extracted_from_count_(
                    geneid2ix, transcript_model, unspliced, cell_bcidx
                )
            if has_onlyintron_model:
                if has_onlyexo_model and not has_mixed_model:
                    return self._extracted_from_count_(
                        geneid2ix, transcript_model, ambiguous, cell_bcidx
                    )
                if not has_onlyexo_model:
                    return self._extracted_from_count_(
                        geneid2ix, transcript_model, ambiguous, cell_bcidx
                    )
            if not has_onlyintron_model and has_onlyexo_model:
                return self._extracted_from_count_(
                    geneid2ix, transcript_model, ambiguous, cell_bcidx
                )
            if has_onlyintron_model:
                return self._extracted_from_count_(
                    geneid2ix, transcript_model, ambiguous, cell_bcidx
                )

    # TODO Rename this here and in `count`
    def _extracted_from_count_(self, geneid2ix, transcript_model, arg2, cell_bcidx):
        # More common situation, normal exonic read, count as spliced
        gene_ix = geneid2ix[transcript_model.geneid]
        arg2[gene_ix, cell_bcidx] += 1
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
    def layers(
        self,
    ) -> list[str]:  # This should be overridden if a different set of layers is desired
        return ["spliced", "unspliced", "ambiguous"]

    @property
    def stranded(self) -> bool:
        return True

    @property
    def perform_validation_markup(self) -> bool:
        return True

    @property
    def accept_discordant(self) -> bool:
        return False

    def count(
        self,
        molitem: Molitem,
        cell_bcidx: int,
        dict_layers_columns: dict[str, np.ndarray],
        geneid2ix: dict[str, int],
    ) -> None:
        if len(molitem.mappings_record) == 0:
            return
        # NOTE This can be simplified qyuite a bit, without loss of acuracy!
        # The hits are not compatible with any annotated transcript model
        spliced = dict_layers_columns["spliced"]
        unspliced = dict_layers_columns["unspliced"]
        ambiguous = dict_layers_columns["ambiguous"]
            # Check that there are not different possible genes ??
        if len({i.geneid for i in molitem.mappings_record.keys()}) == 1:
            gene_check: set[str] = set()

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
                # has_exseg_with_spliced_flag = 0
                has_validated_intron = 0
                has_exin_intron_span = 0
                # has_non3prime = 0
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
                        # if not segment_match.feature.is_last_3prime:
                        #     has_non3prime = 1
                        # if segment_match.is_spliced:
                        #     has_exseg_with_spliced_flag = 1
                if has_validated_intron and not has_exons:
                    has_onlyintron_and_valid_model = 1
                if has_introns and not has_exons:
                    has_onlyintron_model = 1
                if has_exons:
                    if not has_introns:
                        has_onlyexo_model = 1
                    if has_introns:
                        if not has_validated_intron and not has_exin_intron_span:
                            has_invalid_mixed_model = 1
                            has_mixed_model = 1
                        if has_validated_intron and not has_exin_intron_span:
                            has_valid_mixed_model = 1
                            has_mixed_model = 1
                if not has_exin_intron_span:
                    has_only_span_exin_model = 0

            if not multi_gene and not len(molitem.mappings_record) or multi_gene:
                # No gene is compatible with the observation, do not count
                return
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
                # Singleton in validated intron, do not count
                return
            if (
                            has_onlyintron_model
                            and not has_onlyintron_and_valid_model
                            and not has_mixed_model
                            and not has_onlyexo_model
                        ):
                # Singleton in non-validated intron
                return
            if (
                has_invalid_mixed_model
                and not has_valid_mixed_model
                and not has_onlyintron_model
                and not has_onlyexo_model
            ):
                # Not validated and mapping to exon and introns, happens rarely in 10X / irrelevant.
                return
            if (
                has_valid_mixed_model
                and not has_onlyintron_model
                and not has_onlyexo_model
            ):
                # Validated and mapping to exon and introns, happens rarely in 10X. Count as unspliced.
                gene_ix = geneid2ix[transcript_model.geneid]
                unspliced[gene_ix, cell_bcidx] += 1
                return
            if has_onlyintron_model:
                if has_onlyexo_model and not has_mixed_model:
                    # Ambiguity among the transcript models compatible with the mapping, most common case! Count ambiguous
                    gene_ix = geneid2ix[transcript_model.geneid]
                    ambiguous[gene_ix, cell_bcidx] += 1
                    return
                if not has_onlyexo_model:
                    # Ambiguity among the transcript models compatible with the mapping. Very rare.
                    gene_ix = geneid2ix[transcript_model.geneid]
                    ambiguous[gene_ix, cell_bcidx] += 1
                    return
            if not has_onlyintron_model and has_onlyexo_model:
                # Ambiguity among the transcript models compatible with the mapping. Count ambiguous
                gene_ix = geneid2ix[transcript_model.geneid]
                ambiguous[gene_ix, cell_bcidx] += 1
                return
            if has_onlyintron_model:
                # Ambiguity among the transcript models compatible with the mapping. Very rare.
                gene_ix = geneid2ix[transcript_model.geneid]
                ambiguous[gene_ix, cell_bcidx] += 1
                return


class Discordant10X(Logic):
    """Just a test"""

    def __init__(self) -> None:
        self.name = "Discordant10X"

    @property
    def layers(
        self,
    ) -> list[str]:  # This should be overridden if a different set of layers is desired
        return ["spliced", "unspliced", "ambiguous"]

    @property
    def stranded(self) -> bool:
        return True

    @property
    def perform_validation_markup(self) -> bool:
        return True

    @property
    def accept_discordant(self) -> bool:
        return True

    def count(
        self,
        molitem: Molitem,
        cell_bcidx: int,
        dict_layers_columns: dict[str, np.ndarray],
        geneid2ix: dict[str, int],
    ) -> None:
        if len(molitem.mappings_record) == 0:
            return
        # NOTE This can be simplified qyuite a bit, without loss of acuracy!
        # The hits are not compatible with any annotated transcript model
        spliced = dict_layers_columns["spliced"]
        unspliced = dict_layers_columns["unspliced"]
        ambiguous = dict_layers_columns["ambiguous"]
            # Check that there are not different possible genes ??
        if len({i.geneid for i in molitem.mappings_record.keys()}) == 1:
            gene_check: set[str] = set()

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
                # has_exseg_with_spliced_flag = 0
                has_validated_intron = 0
                has_exin_intron_span = 0
                # has_non3prime = 0
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
                        # if not segment_match.feature.is_last_3prime:
                        #     has_non3prime = 1
                        # if segment_match.is_spliced:
                        #     has_exseg_with_spliced_flag = 1
                if has_validated_intron and not has_exons:
                    has_onlyintron_and_valid_model = 1
                if has_introns and not has_exons:
                    has_onlyintron_model = 1
                if has_exons:
                    if not has_introns:
                        has_onlyexo_model = 1
                    if has_introns:
                        if not has_validated_intron and not has_exin_intron_span:
                            has_invalid_mixed_model = 1
                            has_mixed_model = 1
                        if has_validated_intron and not has_exin_intron_span:
                            has_valid_mixed_model = 1
                            has_mixed_model = 1
                if not has_exin_intron_span:
                    has_only_span_exin_model = 0

            if not multi_gene and not len(molitem.mappings_record) or multi_gene:
                # No gene is compatible with the observation, do not count
                return
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
                # Singleton in validated intron
                gene_ix = geneid2ix[transcript_model.geneid]
                unspliced[gene_ix, cell_bcidx] += 1
                return
            if (
                            has_onlyintron_model
                            and not has_onlyintron_and_valid_model
                            and not has_mixed_model
                            and not has_onlyexo_model
                        ):
                # Singleton in non-validated intron
                gene_ix = geneid2ix[transcript_model.geneid]
                unspliced[gene_ix, cell_bcidx] += 1
                return
            if (
                has_invalid_mixed_model
                and not has_valid_mixed_model
                and not has_onlyintron_model
                and not has_onlyexo_model
            ):
                # Not validated and mapping to exon and introns, happens rarely in 10X / irrelevant. Count anyways
                gene_ix = geneid2ix[transcript_model.geneid]
                unspliced[gene_ix, cell_bcidx] += 1
                return
            if (
                has_valid_mixed_model
                and not has_onlyintron_model
                and not has_onlyexo_model
            ):
                # Validated and mapping to exon and introns, happens rarely in 10X. Count as unspliced.
                gene_ix = geneid2ix[transcript_model.geneid]
                unspliced[gene_ix, cell_bcidx] += 1
                return
            if has_onlyintron_model:
                if has_onlyexo_model and not has_mixed_model:
                    # Ambiguity among the transcript models compatible with the mapping, most common case! Count ambiguous
                    gene_ix = geneid2ix[transcript_model.geneid]
                    ambiguous[gene_ix, cell_bcidx] += 1
                    return
                if not has_onlyexo_model:
                    # Ambiguity among the transcript models compatible with the mapping. Very rare. Count ambiguous
                    gene_ix = geneid2ix[transcript_model.geneid]
                    ambiguous[gene_ix, cell_bcidx] += 1
                    return
            if not has_onlyintron_model and has_onlyexo_model:
                # Ambiguity among the transcript models compatible with the mapping. Count ambiguous
                gene_ix = geneid2ix[transcript_model.geneid]
                ambiguous[gene_ix, cell_bcidx] += 1
                return
            if has_onlyintron_model:
                # Ambiguity among the transcript models compatible with the mapping. Very rare. Count ambiguous
                gene_ix = geneid2ix[transcript_model.geneid]
                ambiguous[gene_ix, cell_bcidx] += 1
                return


class SmartSeq2(Logic):
    """SmartSeq2 logic"""

    def __init__(self) -> None:
        self.name = "SmartSeq2"

    @property
    def layers(
        self,
    ) -> list[str]:  # This should be overridden if a different set of layers is desired
        return ["spliced", "unspliced", "ambiguous", "spanning"]

    @property
    def stranded(self) -> bool:
        return False

    @property
    def perform_validation_markup(self) -> bool:
        return False

    @property
    def accept_discordant(self) -> bool:
        return False

    def count(
        self,
        molitem: Molitem,
        cell_bcidx: int,
        dict_layers_columns: dict[str, np.ndarray],
        geneid2ix: dict[str, int],
    ) -> None:
        if len(molitem.mappings_record) == 0:
            return
        # NOTE This can be simplified qyuite a bit, without loss of acuracy!
        # The hits are not compatible with any annotated transcript model
        spliced = dict_layers_columns["spliced"]
        unspliced = dict_layers_columns["unspliced"]
        ambiguous = dict_layers_columns["ambiguous"]
        spanning = dict_layers_columns["spanning"]

            # Check that there are not different possible genes ??
        if len({i.geneid for i in molitem.mappings_record.keys()}) == 1:
            gene_check: set[str] = set()

            has_onlyintron_model = 0
            has_only_span_exin_model = 1
            has_onlyexo_model = 0
            has_mixed_model = 0
            multi_gene = 0
            for transcript_model, segments_list in molitem.mappings_record.items():
                gene_check.add(transcript_model.geneid)
                if len(gene_check) > 1:
                    multi_gene = 1
                has_introns = 0
                has_exons = 0
                # has_exseg_with_spliced_flag = 0
                has_exin_intron_span = 0
                for segment_match in segments_list:
                    if segment_match.maps_to_intron:
                        has_introns = 1
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
                        # if segment_match.is_spliced:
                        #     has_exseg_with_spliced_flag = 1
                if has_introns and not has_exons:
                    has_onlyintron_model = 1
                if has_exons:
                    if not has_introns:
                        has_onlyexo_model = 1
                    if has_introns and not has_exin_intron_span:
                        # has_valid_mixed_model = 1
                        has_mixed_model = 1
                if not has_exin_intron_span:
                    has_only_span_exin_model = 0

            if not multi_gene and not len(molitem.mappings_record) or multi_gene:
                # NOTE it does not happen for Smartseq2
                # No gene is compatible with the observation, do not count
                return
            if has_onlyexo_model and not has_onlyintron_model and not has_mixed_model:
                # More common situation, normal exonic read, count as spliced
                gene_ix = geneid2ix[transcript_model.geneid]
                spliced[gene_ix, cell_bcidx] += 1
                return
            if has_only_span_exin_model:
                # NOTE This is what I want to count as spanning
                # All the compatible transcript models have spanning exon-intron boundaries, count unspliced
                gene_ix = geneid2ix[transcript_model.geneid]
                spanning[gene_ix, cell_bcidx] += 1
                return
            if has_onlyintron_model and not has_mixed_model:
                if not has_onlyexo_model:
                    gene_ix = geneid2ix[transcript_model.geneid]
                    unspliced[gene_ix, cell_bcidx] += 1
                    return
                # Ambiguity among the transcript models compatible with the mapping, most common case! Count ambiguous
                gene_ix = geneid2ix[transcript_model.geneid]
                ambiguous[gene_ix, cell_bcidx] += 1
                return
            if not has_onlyintron_model and has_onlyexo_model:
                # NOTE has_mixed model is used only here in this logic
                # Ambiguity among the transcript models compatible with the mapping. Count ambiguous
                gene_ix = geneid2ix[transcript_model.geneid]
                ambiguous[gene_ix, cell_bcidx] += 1
                return


Default: type = Permissive10X
