from pathlib import Path
from typing import Iterator, Optional, Union

class NeedletailError(Exception):
    """
    Raised when parsing invalid fastx data with parse_fastx_file or parse_fastx_string
    """

class FastxReader(Iterator[Record]):
    """
    An iterator that yields sequence records.

    Yields
    ------
    Record
        A `Record` object representing a sequence record.

    See also
    --------
    parse_fastx_file:
        A function to parse sequence records from a FASTA/FASTQ file.
    parse_fastx_string:
        A function to parse sequence records from a FASTA/FASTQ string.
    Record:
        A class representing a FASTA/FASTQ sequence record.
    """

class Record:
    """
    A record representing a biological sequence.

    Parameters
    ----------
    id : str
        The identifier of the sequence record.
    seq : str
        A string representing the sequence.

    Attributes
    ----------
    id : str
        The identifier of the sequence record. In a FASTA file, this is the
        string containing all characters (including whitespaces) after the
        leading '>' character. In a FASTQ file, this is the string containing
        all characters (including whitespaces) after the leading '@' character.
    seq : str
        A string representing the sequence.
    qual : str, optional
        A string representing the quality scores of the sequence. If the object
        represents a FASTA record, this attribute will be `None`.
    name : str
        The name of the sequence record. This is the string before the first
        whitespace character in the `id` attribute.
    description : str, optional
        The description of the sequence record. This is the string after the
        first whitespace character in the `id` attribute. If the `id` attribute
        contains no whitespace characters, this attribute will be `None`.

    Methods
    -------
    is_fasta
        Check if the object represents a FASTA record.
    is_fastq
        Check if the object represents a FASTQ record.
    normalize(iupac)
        Normalize the sequence stored in the `seq` attribute of the object.
    kmers(k)
        Return an iterator over the k-mers of the sequence stored in the `seq`
        attribute.
    canonical_kmers(k)
        Return an iterator over the canonical k-mers of the sequence stored in
        the `seq` attribute.
    bit_kmers(k, canonical)
        Return an iterator over the k-mers of the sequence stored in the `seq`
        attribute as 2-bit packed integers.
    """

    id: str
    seq: str
    name: str
    description: Optional[str]
    qual: Optional[str]

    def is_fasta(self) -> bool:
        """
        Check if the object represents a FASTA record.

        Returns
        -------
        bool
            `True` if the record lacks quality information, otherwise `False`.
        """
        pass

    def is_fastq(self) -> bool:
        """
        Check if the object represents a FASTQ record.

        Returns
        -------
        bool
            `True` if the record has quality information, otherwise `False`.
        """
        pass

    def normalize(self, iupac: bool) -> None:
        """
        Normalize the sequence stored in the `seq` attribute of the object.

        See also
        --------
        normalize_seq: A function to normalize nucleotide sequence strings.

        Notes
        -----
        The `normalize` method is designed for nucleotide sequences only. If
        used with protein sequences, it will incorrectly process amino acid
        characters as if they were nucleotides.
        """
        pass

    def kmers(self, k: int) -> Iterator[str]:
        """
        Return an iterator over the k-mers of the sequence stored in the `seq`
        attribute.

        See also
        --------
        kmers: A function to split a sequence string into k-mers.
        """
        pass

    def canonical_kmers(self, k: int) -> Iterator[tuple[int, str, bool]]:
        """
        Return an iterator over the canonical k-mers of the sequence stored in
        the `seq` attribute.

        See also
        --------
        canonical_kmers:
            A function to split a nucleotide sequence string into canonical
            k-mers.
        """
        pass

    def bit_kmers(
        self, k: int, canonical: bool = False
    ) -> Iterator[tuple[int, int, bool]]:
        """
        Return an iterator over the k-mers of the sequence stored in the `seq`
        attribute as 2-bit packed integers.

        See also
        --------
        bit_kmers:
            A function to split a nucleotide sequence string into 2-bit packed
            integer k-mers.
        """
        pass

def parse_fastx_file(path: Union[str, Path]) -> FastxReader:
    """
    Returns an iterator that parses a FASTA/FASTQ file and yields sequence
    records.

    Parameters
    ----------
    path : str or pathlib.Path
        The path to a FASTA/FASTQ file.

    Returns
    -------
    FastxReader
        A `FastxReader` iterator that yields `Record` objects representing
        sequences from the input file.

    Raises
    ------
    NeedletailError
        If an error occurs while reading and parsing the input file.

    See also
    --------
    parse_fastx_string:
        A function to parse sequence records from a FASTA/FASTQ string.
    FastxReader:
        A class with instances that are iterators that yield `Record` objects.
    """
    pass

def parse_fastx_string(fastx_string: str) -> FastxReader:
    """
    Returns an iterator that parses a FASTA/FASTQ string and yields sequence
    records.

    Parameters
    ----------
    content : str
        A string containing FASTA/FASTQ-formatted sequence records.

    Returns
    -------
    FastxReader
        A `FastxReader` iterator that yields `Record` objects representing
        sequences from the input string.

    Raises
    ------
    NeedletailError
        If an error occurs while parsing the input string.

    See also
    --------
    parse_fastx_file:
        A function to parse sequence records from a FASTA/FASTQ file.
    FastxReader:
        A class with instances that are iterators that yield `Record` objects.
    """
    pass

def normalize_seq(seq: str, iupac: bool) -> str:
    """
    Normalize the sequence string of nucleotide records by:

    - Converting lowercase characters to uppercase.
    - Removing whitespace and newline characters.
    - Replacing 'U' with 'T'.
    - Replacing '.' and '~' with '-'.
    - Replacing characters not in 'ACGTN-' with 'N', unless `iupac` is `True`,
      in which case characters representing nucleotide ambiguity are not
      replaced.

    Parameters
    ----------
    seq : str
        A string representing a nucleotide sequence.
    iupac : bool, default: False
        If `True`, characters representing nucleotide ambiguity ('B', 'D',
        'H', 'V', 'R', 'Y', 'S', 'W', 'K', and 'M', and their lowercase
        forms) will not be converted to 'N'. Lowercase characters will still
        be converted to uppercase.

    Returns
    -------
    str
        The normalized sequence string.

    Notes
    -----
    The `normalize_seq` function is designed for nucleotide sequences only. If
    used with protein sequences, it will incorrectly process amino acid
    characters as if they were nucleotides.
    """
    pass

def reverse_complement(seq: str) -> str:
    """
    Compute the reverse complement of a nucleotide sequence.

    Parameters
    ----------
    seq : str
        A string representing a nucleotide sequence.

    Returns
    -------
    str
        The reverse complement of the input nucleotide sequence.

    Notes
    -----
    The `reverse_complement` method is designed for nucleotide sequences
    only. If used with protein sequences, it will incorrectly process
    amino acid characters as if they were nucleotides.
    """
    pass

def kmers(seq: str, k: int) -> Iterator[str]:
    """
    Split a sequence into overlapping k-mers.

    Parameters
    ----------
    seq : str
        A string representing a nucleotide or amino acid sequence.
    k : int
        The length of each k-mer (1-255).

    Returns
    -------
    iterator of str
        Every window of length `k` in the sequence, in order. Ambiguous bases
        such as 'N' are not skipped.
    """
    pass

def canonical_kmers(seq: str, k: int) -> Iterator[tuple[int, str, bool]]:
    """
    Split a nucleotide sequence into overlapping canonical k-mers.

    A canonical k-mer is the lexicographically smaller of a k-mer and its
    reverse complement. K-mers containing bases other than 'ACGT' (in either
    case) are skipped.

    Parameters
    ----------
    seq : str
        A string representing a nucleotide sequence.
    k : int
        The length of each k-mer (1-255).

    Returns
    -------
    iterator of tuple of (int, str, bool)
        For each k-mer, its start position in `seq`, the canonical k-mer, and
        whether the canonical k-mer is the reverse complement of the original.

    Notes
    -----
    The `canonical_kmers` function is designed for nucleotide sequences only. The
    input is not validated. With protein sequences, windows consisting only of
    the amino acids 'A', 'C', 'G', and 'T' are treated as nucleotides and the
    rest are skipped, so the result is meaningless rather than an error.
    """
    pass

def bit_kmers(
    seq: str, k: int, canonical: bool = False
) -> Iterator[tuple[int, int, bool]]:
    """
    Split a nucleotide sequence into overlapping k-mers encoded as 2-bit
    packed integers.

    Each base is encoded as two bits (A=0, C=1, G=2, T=3), with the first base
    in the most significant position. K-mers containing bases other than
    'ACGT' (in either case) are skipped.

    Parameters
    ----------
    seq : str
        A string representing a nucleotide sequence.
    k : int
        The length of each k-mer (1-32).
    canonical : bool, default=False
        If `True`, return the smaller of each k-mer and its reverse complement.

    Returns
    -------
    iterator of tuple of (int, int, bool)
        For each k-mer, its start position in `seq`, the packed integer, and
        whether the returned k-mer is the reverse complement of the original
        (always `False` unless `canonical` is `True`).

    Notes
    -----
    The `bit_kmers` function is designed for nucleotide sequences only. The
    input is not validated. With protein sequences, windows consisting only of
    the amino acids 'A', 'C', 'G', and 'T' are treated as nucleotides and the
    rest are skipped, so the result is meaningless rather than an error.
    """
    pass

def decode_phred(qual: str, base_64: bool) -> tuple[int]:
    """
    Decode Phred quality strings to quality scores.

    Parameters
    ----------
    phred : str
        A string representing Phred-encoded quality strings.
    base_64 : bool, default=False
        If `True`, return the quality using the Phred+64 encoding, otherwise
        the Phred+33 encoding will be used.

    Returns
    -------
    tuple of int
        A list of integers representing quality scores derived from the
        probability of a base-calling error using a logarithmic transformation.
    """
    pass
