use num_integer::Roots;

#[derive(Debug)]
pub enum PairingError {
    Overflow,
}

pub fn compute_cantor_pair(l: u64, r: u64) -> Result<u64, PairingError> {
    let sum = l.checked_add(r).ok_or(PairingError::Overflow)?;
    let product = sum.checked_mul(sum + 1).ok_or(PairingError::Overflow)?;
    let half = product.checked_div(2).ok_or(PairingError::Overflow)?;
    half.checked_add(r).ok_or(PairingError::Overflow)
}

pub fn compute_cantor_unpair(p: u64) -> Result<(u64, u64), PairingError> {
    let w = ((8 * p + 1).sqrt() - 1) / 2;
    let t = (w * (w + 1)) / 2;
    let r = p - t;
    let l = w - r;
    Ok((l, r))
}

pub fn compute_peter_pair(l: u64, r: u64) -> Result<u64, PairingError> {
    let (shell, step) = if l > r { (l, r) } else { (r, l) };
    let flag = if step == r { 0 } else { 1 };
    let shell_squared = shell.checked_mul(shell).ok_or(PairingError::Overflow)?;
    let step_doubled = step.checked_mul(2).ok_or(PairingError::Overflow)?;
    shell_squared
        .checked_add(step_doubled)
        .and_then(|v| v.checked_add(flag))
        .ok_or(PairingError::Overflow)
}

pub fn compute_peter_unpair(p: u64) -> Result<(u64, u64), PairingError> {
    let shell = p.sqrt();
    let step = (p - shell.pow(2)) / 2;
    if p % 2 == 0 {
        Ok((shell, step))
    } else {
        Ok((step, shell))
    }
}

pub fn compute_rosenberg_strong_pair(l: u64, r: u64) -> Result<u64, PairingError> {
    let shell = l.max(r);
    let shell_squared = shell.checked_mul(shell).ok_or(PairingError::Overflow)?;
    shell_squared
        .checked_add(shell)
        .and_then(|v| v.checked_add(l))
        .and_then(|v| v.checked_sub(r))
        .ok_or(PairingError::Overflow)
}

pub fn compute_rosenberg_strong_unpair(p: u64) -> Result<(u64, u64), PairingError> {
    let shell = p.sqrt();
    let step = p - shell.pow(2);
    if step < shell {
        Ok((step, shell))
    } else {
        Ok((shell, 2 * shell - step))
    }
}

pub fn compute_szudzik_pair(l: u64, r: u64) -> Result<u64, PairingError> {
    if l != l.max(r) {
        let r_squared = r.checked_mul(r).ok_or(PairingError::Overflow)?;
        r_squared.checked_add(l).ok_or(PairingError::Overflow)
    } else {
        let l_squared = l.checked_mul(l).ok_or(PairingError::Overflow)?;
        l_squared
            .checked_add(l)
            .and_then(|v| v.checked_add(r))
            .ok_or(PairingError::Overflow)
    }
}

pub fn compute_szudzik_unpair(p: u64) -> Result<(u64, u64), PairingError> {
    let sqrt_p = p.sqrt();
    if p - sqrt_p.pow(2) < sqrt_p {
        Ok((p - sqrt_p.pow(2), sqrt_p))
    } else {
        Ok((sqrt_p, p - sqrt_p.pow(2) - sqrt_p))
    }
}

pub fn compute_hagen_pair(l: u64, r: u64) -> Result<u64, PairingError> {
    let (shell, step) = if l > r { (l, r) } else { (r, l) };
    let flag = match (shell % 2, step) {
        (0, s) if s == l => 0,
        (1, s) if s == r => 0,
        _ => 1,
    };
    let shell_squared = shell.checked_mul(shell).ok_or(PairingError::Overflow)?;
    let step_doubled = step.checked_mul(2).ok_or(PairingError::Overflow)?;
    shell_squared
        .checked_add(step_doubled)
        .and_then(|v| v.checked_add(flag))
        .ok_or(PairingError::Overflow)
}

pub fn compute_hagen_unpair(p: u64) -> Result<(u64, u64), PairingError> {
    let shell = p.sqrt();
    let step = (p - shell.pow(2)) / 2;
    if p % 2 == 0 {
        Ok((step, shell))
    } else {
        Ok((shell, step))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cantor_pair() {
        let pairs = vec![
            (1958, 1962),
            (1958, 1970),
            (1958, 1994),
            (1958, 2002),
            (1962, 1970),
            (1962, 1994),
            (1962, 2002),
            (1970, 1994),
            (1970, 2002),
            (1994, 2002),
        ];
        let expected = vec![
            7687122, 7718526, 7813122, 7844782, 7734248, 7828940, 7860632, 7860624, 7892380,
            7988008,
        ];
        for ((l, r), expected) in pairs.iter().zip(expected.iter()) {
            let p = compute_cantor_pair(*l, *r).unwrap();
            assert_eq!(p, *expected);
        }
    }

    #[test]
    fn test_cantor_unpair() {
        let pairs = vec![
            (1958, 1962),
            (1958, 1970),
            (1958, 1994),
            (1958, 2002),
            (1962, 1970),
            (1962, 1994),
            (1962, 2002),
            (1970, 1994),
            (1970, 2002),
            (1994, 2002),
        ];
        let expected = vec![
            7687122, 7718526, 7813122, 7844782, 7734248, 7828940, 7860632, 7860624, 7892380,
            7988008,
        ];
        for ((l, r), expected) in pairs.iter().zip(expected.iter()) {
            let (ul, ur) = compute_cantor_unpair(*expected).unwrap();
            assert_eq!((*l, *r), (ul, ur));
        }
    }

    #[test]
    fn test_peter_pair() {
        let pairs = vec![
            (1958, 1962),
            (1958, 1970),
            (1958, 1994),
            (1958, 2002),
            (1962, 1970),
            (1962, 1994),
            (1962, 2002),
            (1970, 1994),
            (1970, 2002),
            (1994, 2002),
        ];
        let expected = vec![
            3853361, 3884817, 3979953, 4011921, 3884825, 3979961, 4011929, 3979977, 4011945,
            4011993,
        ];
        for ((l, r), expected) in pairs.iter().zip(expected.iter()) {
            let p = compute_peter_pair(*l, *r).unwrap();
            assert_eq!(p, *expected);
        }
    }

    #[test]
    fn test_peter_unpair() {
        let pairs = vec![
            (1958, 1962),
            (1958, 1970),
            (1958, 1994),
            (1958, 2002),
            (1962, 1970),
            (1962, 1994),
            (1962, 2002),
            (1970, 1994),
            (1970, 2002),
            (1994, 2002),
        ];
        let expected = vec![
            3853361, 3884817, 3979953, 4011921, 3884825, 3979961, 4011929, 3979977, 4011945,
            4011993,
        ];
        for ((l, r), expected) in pairs.iter().zip(expected.iter()) {
            let (ul, ur) = compute_peter_unpair(*expected).unwrap();
            assert_eq!((*l, *r), (ul, ur));
        }
    }

    #[test]
    fn test_rosenberg_strong_pair() {
        let pairs = vec![
            (1958, 1962),
            (1958, 1970),
            (1958, 1994),
            (1958, 2002),
            (1962, 1970),
            (1962, 1994),
            (1962, 2002),
            (1970, 1994),
            (1970, 2002),
            (1994, 2002),
        ];
        let expected = vec![
            3851402, 3882858, 3977994, 4009962, 3882862, 3977998, 4009966, 3978006, 4009974,
            4009998,
        ];
        for ((l, r), expected) in pairs.iter().zip(expected.iter()) {
            let p = compute_rosenberg_strong_pair(*l, *r).unwrap();
            assert_eq!(p, *expected);
        }
    }

    #[test]
    fn test_rosenberg_strong_unpair() {
        let pairs = vec![
            (1958, 1962),
            (1958, 1970),
            (1958, 1994),
            (1958, 2002),
            (1962, 1970),
            (1962, 1994),
            (1962, 2002),
            (1970, 1994),
            (1970, 2002),
            (1994, 2002),
        ];
        let expected = vec![
            3851402, 3882858, 3977994, 4009962, 3882862, 3977998, 4009966, 3978006, 4009974,
            4009998,
        ];
        for ((l, r), expected) in pairs.iter().zip(expected.iter()) {
            let (ul, ur) = compute_rosenberg_strong_unpair(*expected).unwrap();
            assert_eq!((*l, *r), (ul, ur));
        }
    }

    #[test]
    fn test_szudzik_pair() {
        let pairs = vec![
            (1958, 1962),
            (1958, 1970),
            (1958, 1994),
            (1958, 2002),
            (1962, 1970),
            (1962, 1994),
            (1962, 2002),
            (1970, 1994),
            (1970, 2002),
            (1994, 2002),
        ];
        let expected = vec![
            3851402, 3882858, 3977994, 4009962, 3882862, 3977998, 4009966, 3978006, 4009974,
            4009998,
        ];
        for ((l, r), expected) in pairs.iter().zip(expected.iter()) {
            let p = compute_szudzik_pair(*l, *r).unwrap();
            assert_eq!(p, *expected);
        }
    }

    #[test]
    fn test_szudzik_unpair() {
        let pairs = vec![
            (1958, 1962),
            (1958, 1970),
            (1958, 1994),
            (1958, 2002),
            (1962, 1970),
            (1962, 1994),
            (1962, 2002),
            (1970, 1994),
            (1970, 2002),
            (1994, 2002),
        ];
        let expected = vec![
            3851402, 3882858, 3977994, 4009962, 3882862, 3977998, 4009966, 3978006, 4009974,
            4009998,
        ];
        for ((l, r), expected) in pairs.iter().zip(expected.iter()) {
            let (ul, ur) = compute_szudzik_unpair(*expected).unwrap();
            assert_eq!((*l, *r), (ul, ur));
        }
    }

    #[test]
    fn test_hagen_pair() {
        let pairs = vec![
            (1958, 1962),
            (1958, 1970),
            (1958, 1994),
            (1958, 2002),
            (1962, 1970),
            (1962, 1994),
            (1962, 2002),
            (1970, 1994),
            (1970, 2002),
            (1994, 2002),
        ];
        let expected = vec![
            3853360, 3884816, 3979952, 4011920, 3884824, 3979960, 4011928, 3979976, 4011944,
            4011992,
        ];
        for ((l, r), expected) in pairs.iter().zip(expected.iter()) {
            let p = compute_hagen_pair(*l, *r).unwrap();
            assert_eq!(p, *expected);
        }
    }

    #[test]
    fn test_hagen_unpair() {
        let pairs = vec![
            (1958, 1962),
            (1958, 1970),
            (1958, 1994),
            (1958, 2002),
            (1962, 1970),
            (1962, 1994),
            (1962, 2002),
            (1970, 1994),
            (1970, 2002),
            (1994, 2002),
        ];
        let expected = vec![
            3853360, 3884816, 3979952, 4011920, 3884824, 3979960, 4011928, 3979976, 4011944,
            4011992,
        ];
        for ((l, r), expected) in pairs.iter().zip(expected.iter()) {
            let (ul, ur) = compute_hagen_unpair(*expected).unwrap();
            assert_eq!((*l, *r), (ul, ur));
        }
    }

    #[test]
    fn test_cantor_overflow() {
        assert!(compute_cantor_pair(u64::MAX, u64::MAX).is_err());
    }

    #[test]
    fn test_peter_overflow() {
        assert!(compute_peter_pair(u64::MAX, u64::MAX).is_err());
    }

    #[test]
    fn test_rosenberg_strong_overflow() {
        assert!(compute_rosenberg_strong_pair(u64::MAX, u64::MAX).is_err());
    }

    #[test]
    fn test_szudzik_overflow() {
        assert!(compute_szudzik_pair(u64::MAX, u64::MAX).is_err());
    }

    #[test]
    fn test_hagen_overflow() {
        assert!(compute_hagen_pair(u64::MAX, u64::MAX).is_err());
    }
}
