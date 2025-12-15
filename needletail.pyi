from pathlib import Path
from typing import Iterator, Optional, Union

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
