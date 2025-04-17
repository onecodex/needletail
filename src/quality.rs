use crate::errors::PhredOffsetError;

/// Encoding for quality scores
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum PhredEncoding {
    Phred33,
    Phred64,
}

/// Decodes Phred quality data to quality scores. Each character is decoded
/// by subtracting the offset corresponding to the specified source encoding
/// (PhredEncoding::Phred33 for Phred+33, or PhredEncoding::Phred64 for Phred+64).
/// If the ASCII value of the character is less than the offset, `PhredOffsetError`
/// is returned.
pub fn decode_phred(qual: &[u8], encoding: PhredEncoding) -> Result<Vec<u8>, PhredOffsetError> {
    let offset = match encoding {
        PhredEncoding::Phred33 => b'!',
        PhredEncoding::Phred64 => b'@',
    };
    let mut scores = Vec::with_capacity(qual.len());
    for &q in qual {
        if q < offset {
            return Err(PhredOffsetError { q, offset });
        }
        scores.push(q - offset);
    }
    Ok(scores)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_decode_phred33() {
        assert_eq!(
            decode_phred(b"#</</BBFFFBF<", PhredEncoding::Phred33),
            Ok(vec![2, 27, 14, 27, 14, 33, 33, 37, 37, 37, 33, 37, 27])
        );
    }

    #[test]
    fn test_decode_phred64() {
        assert_eq!(
            decode_phred(b"B[N[Naaeeeae[", PhredEncoding::Phred64),
            Ok(vec![2, 27, 14, 27, 14, 33, 33, 37, 37, 37, 33, 37, 27])
        );
    }

    #[test]
    fn test_decode_phred33_error() {
        assert_eq!(
            decode_phred(b"#</</BBFFFBF ", PhredEncoding::Phred33),
            Err(PhredOffsetError { q: 32, offset: 33 })
        );
    }

    #[test]
    fn test_decode_phred64_error() {
        assert_eq!(
            decode_phred(b"B[N[Naaeeeae?", PhredEncoding::Phred64),
            Err(PhredOffsetError { q: 63, offset: 64 })
        );
    }
}
