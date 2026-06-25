import re
import time
from collections import Counter
from collections import defaultdict
from fractions import Fraction
import itertools

start_time = time.perf_counter()

commutator = '[((s1-1)*Y1*Z1 + (s2-1)*Y2*Z2 + (s3-1)*Y3*Z3)*(Z1*Y1 + Z2*Y2 + Z3*Y3)*(Z1^n1)*(Z2^n2)*(Z3^n3), a1^m*P1_m]'
#commutator = '[(Y3)*(Z1^n1)*(Z2^n2)*(Z3^(n3+1))*((s2+s3-2)*Y1*Z1 + (s3+s1-2)*Y2*Z2 + (s3-1)*Y3*Z3), a1^m*P1_m]'
#commutator = '[(Y1)*(Y2)*(Y3)*(Z1^p1)*(Z2^p2)*(Z3^p3), a1^m*P1_m]'
#commutator = '[(Y1)*(Y1)*(Y1)*(Y1)*(Y2)*(Y2)*(Y2)*(Y2)*(Y3)*(Y3)*(Y3)*(Y3)*(Z1^p1)*(Z2^p2)*(Z3^p3), a1^m*P1_m]'

z_term_comm = '[(1/l^2)*(A)*(Z1^(n1+1))*(Z2^(n2+1))*(Z3^(n3+1)), a1^m*P1_m]'
#z_term_comm = '[(A*Y1*Z1^(p1+1)*Z2^p2*Z3^p3 + B*Y2*Z1^p1*Z2^(p2+1)*Z3^p3 + C*Y3*Z1^p1*Z2^p2*Z3^(p3+1)), a1^m*P1_m]'

def top_level_split(expression, delimiter='+'):
    """Split the expression by the delimiter only at top-level (depth=0)."""
    parts = []
    bracket_level = 0
    current = []
    for ch in expression:
        if ch == '(':
            bracket_level += 1
            current.append(ch)
        elif ch == ')':
            bracket_level -= 1
            current.append(ch)
        elif ch == delimiter and bracket_level == 0:
            # We found a top-level delimiter
            part = ''.join(current).strip()
            if part:
                parts.append(part)
            current = []
        else:
            current.append(ch)
    # Add the last part
    last_part = ''.join(current).strip()
    if last_part:
        parts.append(last_part)
    return parts

def parse_exponent(exp):
    """
    Same signature as before. Parses the exponent string, returning:
      parenthesized (bool), symbol (str), int_part (int or None)
    """
    exp = exp.strip()
    parenthesized = exp.startswith('(') and exp.endswith(')')
    if parenthesized:
        exp = exp[1:-1].strip()

    # Try matching a trailing +digits
    pattern = re.compile(r'^([A-Za-z0-9]+)(?:\+(\d+))?$')
    m = pattern.match(exp)
    if m:
        symbol = m.group(1)
        int_part = m.group(2)
        if int_part is not None:
            int_part = int(int_part)
        return parenthesized, symbol, int_part
    else:
        # Try matching a trailing -digits
        pattern_minus = re.compile(r'^([A-Za-z0-9]+)([-]\d+)$')
        m2 = pattern_minus.match(exp)
        if m2:
            symbol = m2.group(1)
            int_part = int(m2.group(2))
            return parenthesized, symbol, int_part

        # Otherwise no integer part
        return parenthesized, exp, None

def build_exponent(parenthesized, symbol, int_part):
    """
    Rebuilds exponent from (parenthesized, symbol, int_part), ensuring:
      - no leading "+" for a numeric exponent
      - parentheses used only if there's an internal +/- or originally parenthesized
      - no extraneous spaces
    
    Examples:
      symbol='n3', int_part=2  -> "n3+2"
      symbol='n3', int_part=-1 -> "n3-1"
      symbol='', int_part=3    -> "3"   (omitting leading '+')
      symbol='', int_part=-2   -> "-2"
      symbol='', int_part=0    -> ""    (meaning exponent=0 => effectively 1)
    """

    # 1) Build the "core exponent" string
    if int_part is not None:
        # If there's no symbol, we'll produce just the numeric part (avoid leading '+')
        if symbol == '':
            # purely numeric exponent
            if int_part > 0:
                exp_str = str(int_part)   # e.g. "3", not "+3"
            elif int_part == 0:
                # exponent=0 => effectively "0", but usually means skip the factor
                # We'll return "" here, so the caller can interpret "Z^( )" => skip
                exp_str = ''
            else:
                # negative => e.g. "-2"
                exp_str = str(int_part)
        else:
            # we do have a symbol
            if int_part > 0:
                # e.g. "n3+2"
                exp_str = f"{symbol}+{int_part}"
            elif int_part == 0:
                # effectively just the symbol
                exp_str = symbol
            else:
                # negative offset => e.g. "n3-1"
                exp_str = f"{symbol}{int_part}"
    else:
        # no int_part => just the symbol
        exp_str = symbol

    # 2) Check if exponent is empty => means exponent=0 => return ""
    #    That tells the caller "there's no exponent" => factor is base^(0)->1 => skip
    if exp_str == '':
        return ''

    # 3) Decide if we need parentheses
    #    We only wrap if there's a plus/minus beyond the first character,
    #    or if it was originally parenthesized. 
    #    e.g. "n3+2" => "(n3+2)"
    #         "3"    => no parentheses
    #         "-2"   => first char is '-', so no internal sign => no parentheses
    #         "n3-2" => there's a minus in the 1st char after 'n3'? Actually we check the final string
    # We'll skip the first character in the search for plus/minus
    if len(exp_str) > 1:
        has_inner_sign = any(ch in '+-' for ch in exp_str[1:])
    else:
        has_inner_sign = False

    if has_inner_sign:
        exp_str = f"({exp_str})"
    elif parenthesized:
        exp_str = f"({exp_str})"

    return exp_str

def adjust_exponent(exp, delta):
    """
    Same signature. Adds 'delta' to the exponent's integer part if present,
    then rebuilds via build_exponent.
    """

    parenthesized, symbol, int_part = parse_exponent(exp)
    if int_part is not None:
        new_int = int_part + delta
        if new_int == 0:
            int_part = None
        else:
            int_part = new_int
    else:
        # no int_part => attach '+delta' or '-delta'
        if delta > 0:
            symbol = f"{symbol}+{delta}"
        elif delta < 0:
            symbol = f"{symbol}{delta}"

    if int_part == None and parenthesized:
        parenthesized = False

    return build_exponent(parenthesized, symbol, int_part)

def decrement_exponent(exp):
    return adjust_exponent(exp, -1)

def reorder_factors_in_equation(equation_string):
    """
    Takes an equation string, e.g.:
      '-(1/l^2)*2*Z1^(n1+1)*Z2^(n2+1)*Z3^(n3)*Y2'
    and reorders each term's factors so that:
      1) pure numeric digits (like '2','-3','3.5') come first,
      2) strictly numeric parenthesized fractions (like '(3/2)','(-1.5/4)') next,
      3) leftover factors (non-numeric, non-s/z/y) next,
      4) s-variables (s\d+) next,
      5) Z-variables (Z\d+) next,
      6) Y-variables (Y\d+) last.
    No exponent merging or changes, just reordering.

    Steps:
      - top_level_split(equation_string) => list of terms (with +/- sign).
      - For each term:
          * parse sign,
          * split by '*',
          * reorder factors by the above logic,
          * rejoin with '*',
          * reattach sign if '-'.
      - Rejoin terms with ' + ', skipping a leading '+' on the first if it’s positive.
    """

    # 1) We assume you have a function top_level_split(equation_string)
    terms = top_level_split(equation_string)

    def reorder_factors_in_term(term_body):
        """
        Splits 'term_body' by '*', reorders factors according to:
          group= (priority..., alpha-strings)
           (0, 0) => purely digit numeric e.g. '2','-3','3.5'
           (0, 1) => strictly numeric parenthesized fraction e.g. '(3/2)','(-1.5/4.2)'
           (1,)   => leftover (non-numeric, non-s/z/y)
           (2,)   => s\d+ 
           (3,)   => Z\d+
           (4,)   => Y\d+
        Then sorts by that group, tie => alphabetical.
        """

        raw_factors = [f.strip() for f in term_body.split('*') if f.strip()]

        # 1) purely_digit_pattern => ^[+\-]?\d+(\.\d+)?$
        purely_digit_pattern = re.compile(r'^[+\-]?\d+(?:\.\d+)?$')

        # 2) fraction_pattern => ^\(\s*[+\-]?\d+(\.\d+)?\s*/\s*[+\-]?\d+(\.\d+)?\s*\)$
        #    ensures we only match something like "(3/2)", "(1.5/4.7)", "-(2/3)", etc.
        fraction_pattern = re.compile(r"""
            ^               # start
            \(\s*[+\-]?\d+(?:\.\d+)?    # optional sign, digits, optional .digits
            \s*/\s*
            [+\-]?\d+(?:\.\d+)?         # optional sign, digits, optional .digits
            \s*\)$
            """, re.VERBOSE)

        def factor_key(ff):
            s = ff.strip()
            # group=0 => numeric, sub-group=0 => purely digit
            if purely_digit_pattern.match(s):
                return (0, 0, s)
            # group=0 => numeric, sub-group=1 => fraction
            if fraction_pattern.match(s):
                return (0, 1, s)

            # group=2 => s\d+
            if s.startswith('s'):
                return (2, s)
            # group=3 => Z\d+
            if s.startswith('Z'):
                return (3, s)
            # group=4 => Y\d+
            if s.startswith('Y'):
                return (4, s)

            # leftover => group=1
            return (1, s)

        # sort raw_factors by that key
        raw_factors.sort(key=factor_key)
        return '*'.join(raw_factors)

    reordered_terms = []
    for t in terms:
        s = t.strip()
        sign = ''
        if s.startswith('-'):
            sign='-'
            s=s[1:].strip()
        elif s.startswith('+'):
            sign='+'
            s=s[1:].strip()

        # reorder factors
        reordered = reorder_factors_in_term(s)

        # reattach sign if needed
        if sign=='-':
            reordered_terms.append('-'+reordered)
        elif sign=='+':
            reordered_terms.append('+'+reordered)
        else:
            reordered_terms.append(reordered)

    # rejoin with ' + '
    final_list = []
    for i,rt in enumerate(reordered_terms):
        st = rt.strip()
        if i==0 and st.startswith('+'):
            st = st[1:].strip()
        final_list.append(st)

    return ' + '.join(final_list)

def combine_powers_in_terms(equation_string):
    """
    1) Splits the equation string (assumed to be a sum of terms) using top_level_split.
    2) For each term:
       - Removes any leading '+' or '-' sign.
       - Splits by '*' into factors.
       - For factors whose base is one of: Z.., Y.., s.., or m.., merges repeated occurrences by summing their exponents.
       - Reorders factors so that numeric factors come first, then non-combinable ones,
         then m–factors, s–factors, Z–factors, and finally Y–factors.
       - In the output, if an exponent is 1 it is omitted, if 0 the factor is skipped,
         and for a symbolic exponent exactly one set of parentheses is used.
    3) Reattaches the original sign to each term and joins the terms using ' + '.
    4) Calls reorder_factors_in_equation for any final adjustments.
    """
    # 1) Split into terms
    terms = top_level_split(equation_string)

    ###############################################
    # A) parse_factor: return (base, exponent_str) if factor is of the form X\d+^(...),
    #    otherwise (factor, None)
    ###############################################
    def parse_factor(factor):
        factor = factor.strip()
        # Match things like "Z1^(n1+2)" or "Y2^(stuff)" or "m1^(...)" or "s3^(...)" or "n1^(...)".
        mo = re.match(r'^([ZYsmn]\d+)\^(\(?.*?\)?)$', factor)
        if mo:
            return (mo.group(1), mo.group(2))
        # Match a factor without an explicit exponent, e.g. "Z1", "Y2", "m1", etc.
        mo2 = re.match(r'^([ZYsmn]\d+)$', factor)
        if mo2:
            return (mo2.group(1), None)  # exponent=None means exponent=1
        return (factor, None)  # Not combinable

    ###############################################
    # B) unify_exponents: merge two exponent strings if possible.
    ###############################################
    def unify_exponents(old_exp, new_exp):
        """
        Converts None to '1', removes outer parentheses, and tries to
        parse exponents as (symbol, offset). If both exponents have the same symbol
        (or one is purely numeric) it adds their offsets.
        """
        old_str = '1' if old_exp is None else old_exp.strip()
        new_str = '1' if new_exp is None else new_exp.strip()

        def strip_parens(s):
            s = s.strip()
            if s.startswith('(') and s.endswith(')'):
                return s[1:-1].strip()
            return s

        old_str = strip_parens(old_str)
        new_str = strip_parens(new_str)

        # Try to parse an exponent in the form: symbol+offset (or just a number)
        def parse_simple(e):
            if re.match(r'^[+\-]?\d+$', e):
                return ('', int(e))
            plus_pat = re.match(r'^([A-Za-z0-9]+)\+(\d+)$', e)
            if plus_pat:
                return (plus_pat.group(1), int(plus_pat.group(2)))
            minus_pat = re.match(r'^([A-Za-z0-9]+)(-\d+)$', e)
            if minus_pat:
                return (minus_pat.group(1), int(minus_pat.group(2)))
            return (e, 0)

        def build_simple(sym, off):
            if sym == '' and off == 0:
                return '0'
            elif sym == '' and off != 0:
                return str(off)
            elif off == 0:
                return sym
            elif off > 0:
                return f"{sym}+{off}"
            else:
                return f"{sym}{off}"

        s1, o1 = parse_simple(old_str)
        s2, o2 = parse_simple(new_str)

        if s1 == s2:
            new_sym = s1
        elif s1 == '':
            new_sym = s2
        elif s2 == '':
            new_sym = s1
        else:
            # If the symbols differ, we cannot combine; return the old exponent.
            return old_exp

        merged_off = (o1 or 0) + (o2 or 0)
        return build_simple(new_sym, merged_off)

    ###############################################
    # C) combine_factors_in_term: merges repeated factors in a single term.
    ###############################################
    def combine_factors_in_term(term_body):
        """
        Splits the term into factors by '*', merges repeated occurrences of combinable
        factors (whose base is one of s, Z, Y, or m), then reorders the factors.
        """
        raw_factors = [x.strip() for x in term_body.split('*') if x.strip()]

        store = {}  # key: base, value: merged exponent string
        non_combinable = []

        for fac in raw_factors:
            base, exp_str = parse_factor(fac)
            # For combinable factors (s*, Z*, Y*, or m*), merge their exponents.
            if base.startswith('s') or base.startswith('Z') or base.startswith('Y') or base.startswith('m'):
                if base not in store:
                    store[base] = exp_str
                else:
                    store[base] = unify_exponents(store[base], exp_str)
            else:
                non_combinable.append(fac)

        merged_list = []
        for base, e_str in store.items():
            if e_str is None:
                merged_list.append(base)
                continue
            s = e_str.strip()
            if s in ('', '1'):
                merged_list.append(base)
                continue
            if s == '0':
                continue  # skip factor if exponent is 0
            # Remove extra outer parentheses if any, then add one set if symbolic.
            s2 = s
            while s2.startswith('(') and s2.endswith(')'):
                s2 = s2[1:-1].strip()
            if re.match(r'^[+\-]?\d+$', s2):
                merged_list.append(f"{base}^{s2}")
            else:
                merged_list.append(f"{base}^({s2})")

        # Combine merged factors with the non-combinable ones.
        all_factors = merged_list + non_combinable

        # Define a key for ordering factors:
        # Numeric factors come first, then "leftover" factors, then factors in the order:
        # m (group 2), s (group 3), Z (group 4), Y (group 5).
        numeric_pattern = re.compile(r'^[+\-]?\d+(\.\d+)?$|^\(.*?/.*?\)$')
        def factor_key(ff):
            fstr = ff.strip()
            if numeric_pattern.match(fstr):
                return (0, fstr)
            elif fstr.startswith('m'):
                return (2, fstr)
            elif fstr.startswith('s'):
                return (3, fstr)
            elif fstr.startswith('Z'):
                return (4, fstr)
            elif fstr.startswith('Y'):
                return (5, fstr)
            else:
                return (1, fstr)
        
        all_factors.sort(key=factor_key)
        return '*'.join(all_factors)

    # 2) Process each term: remove any leading sign, combine factors, then reattach sign.
    new_terms = []
    for t in terms:
        raw = t.strip()
        sign = ''
        if raw.startswith('-'):
            sign = '-'
            raw = raw[1:].strip()
        elif raw.startswith('+'):
            sign = '+'
            raw = raw[1:].strip()

        merged = combine_factors_in_term(raw)

        if sign:
            new_terms.append(sign + merged)
        else:
            new_terms.append(merged)

    # 3) Rejoin the terms with ' + ', omitting a leading '+' on the first term.
    out_terms = []
    for i, nt in enumerate(new_terms):
        s = nt.strip()
        if i == 0 and s.startswith('+'):
            s = s[1:].strip()
        out_terms.append(s)

    return reorder_factors_in_equation(' + '.join(out_terms))

def pull_minus_signs_to_front(expr):
    """
    Splits expression by " + ".
    For each term:
      - If it starts with "-", we note that it's already negative.
      - Then remove each '*-' from inside, flipping the sign each time.
      - If sign ends up negative, prepend a minus.
    """

    terms = expr.split(' + ')
    new_terms = []

    for term in terms:
        # Detect if the term starts with minus
        # (or possibly multiple minuses, e.g. "--something")
        # We'll consider any leading minus => start_neg = True
        start_neg = bool(re.match(r'^\s*-\s*', term))

        # Start sign as -1 if it began negative, else +1
        sign = -1 if start_neg else 1

        # Strip only the *first* leading minus or plus from the front,
        # if it exists. (If you have multiple minuses, you can adapt.)
        term = re.sub(r'^[+\-]+', '', term, count=1).strip()

        # Now remove any '*-' sequences inside and flip sign each time
        def minus_replacer(_match):
            nonlocal sign
            sign *= -1
            return '*'
        
        new_term = re.sub(r'\*\-', minus_replacer, term)

        # At the end, if sign is negative, put "-" back
        if sign < 0:
            new_term = '-' + new_term
    
        # also check if there is a (-1)* present anywhere
        sub_term_list = new_term.split('*')
        no_minuses = 0
        for var in sub_term_list:
            if var == '(-1)':
                no_minuses += 1
            elif var == '-(-1)':
                continue
        new_sub_term_list = [var for var in sub_term_list if var != '(-1)' and var != '-(-1)']
        if no_minuses % 2 == 1:
            if new_sub_term_list[0][0] == '-':
                new_sub_term_list[0] = new_sub_term_list[0][1:]
            else:
                new_sub_term_list[0] = '-' + new_sub_term_list[0]


        new_terms.append('*'.join(new_sub_term_list))

    return ' + '.join(new_terms)

def find_unused_letter(expr, letters_to_try="bcefghijklmnopqrstuvwxyz"):
    """
    Returns the first letter from 'letters_to_try' that is NOT already used
    in the expression 'expr'. If all letters in 'letters_to_try' are used,
    it raises a ValueError.
    """
    # Find all letters used in the expression (both lowercase and uppercase)
    used_letters = set(re.findall(r'[a-zA-Z]', expr))

    # Pick the first letter not in used_letters
    for letter in letters_to_try:
        if letter not in used_letters:
            return letter
    
    # If we get here, all letters_to_try are in use
    raise ValueError("All candidate letters are already used in the expression.")

def extract_top_level_factors(term, commutator=True):
    if commutator:
        lhs_term = term.split(',')[0][1:]
    else:
        lhs_term = term
    factors = []
    depth = 0
    current = []
    for ch in lhs_term:
        if ch == '(':
            depth += 1
            current.append(ch)
        elif ch == ')':
            depth -= 1
            current.append(ch)
        elif ch == '*' and depth == 0:
            factor = ''.join(current).strip()[1:-1]
            if factor:
                factors.append(factor)
            current = []
        else:
            current.append(ch)
    
    if current:
        last_factor = ''.join(current).strip()[1:-1]
    if last_factor:
        factors.append(last_factor)

    return factors

def expand_factor_list(factors):
    terms_to_expand = []
    for factor in factors:
        terms_in_factor = factor.split(' + ')
        terms_to_expand.append(terms_in_factor)

    while len(terms_to_expand) > 1:
        expand_first_and_second = [terms_to_expand[0][i] + '*' + terms_to_expand[1][j] for i in range(len(terms_to_expand[0])) for j in range(len(terms_to_expand[1]))]
        terms_to_expand = [expand_first_and_second] + terms_to_expand[2:]
    
    return [combine_powers_in_terms(terms_to_expand[0][i]) for i in range(len(terms_to_expand[0]))]

def put_in_commutator(terms):
    # Given a list of terms like:
    # ['(s1-1)*Y1^2*Z1^(n1+2)*Z2^n2*Z3^n3', '(s1-1)*Y1*Z1^(n1+1)*Z2^(n2+1)*Y2*Z3^n3', ...]
    # This function returns:
    # ['[(s1-1)*Y1^2*Z1^(n1+2)*Z2^n2*Z3^n3, a1_m*P1^m]', '[(s1-1)*Y1*Z1^(n1+1)*Z2^(n2+1)*Y2*Z3^n3, a1_m*P1^m]', ...]

    all_terms = ' + '.join(terms)
    new_letter = find_unused_letter(all_terms)

    return [f'[{term}, a1_{new_letter}*P1^{new_letter}]' for term in terms]

def flip_eq_sign(equation):
    terms = top_level_split(equation)
    sign_flipped_terms = []
    for term in terms:
        if term[0] == '-':
            sign_flipped_terms.append(term[1:])
        else:
            sign_flipped_terms.append('-' + term)
        
    return ' + '.join(sign_flipped_terms)


all_comms_to_compute = put_in_commutator(expand_factor_list(extract_top_level_factors(commutator)))
all_z_comms_to_compute = put_in_commutator(expand_factor_list(extract_top_level_factors(z_term_comm)))


def apply_leibniz_rule(comm_expr):
    # comm_expr is something like:
    # "[(s1-1)*Y1^2*Z1^(n1+2)*Z2^n2*Z3^n3, a1_m*P1^m]"
    
    expr = comm_expr[1:-1]

    lhs = expr.split(',')[0]
    rhs = expr.split(',')[1]

    # Split the LHS by '*'
    factors = [f.strip() for f in lhs.split('*') if f.strip()]

    # Identify which factors are variables and which are constants
    # For simplicity, consider constants as those not starting with Y or Z.
    # Variables start with Y or Z.
    # Example:
    # (s1-1)*Y1^2*Z1^(n1+2)*Z2^n2*Z3^n3
    # constants: (s1-1)
    # variables: Y1^2, Z1^(n1+2), Z2^n2, Z3^n3
    constants = []
    variables = []
    for f in factors:
        # Check first character to decide if variable or constant
        if f[0] in ['Y','Z']:
            variables.append(f)
        else:
            constants.append(f)

    # Apply Leibniz rule:
    # [C*V1*V2*...*Vn, X] = sum_{k} C*(V1)*...(V_{k-1})*[V_k, X]*(V_{k+1})*...*Vn
    # where C is product of constants, and V_i are variable factors.
    # constants commute, so factor them out front each time.

    terms = []
    # We'll go from right to left:
    # Actually, the order of Leibniz expansion doesn't matter, we'll produce terms for each variable
    # For variable at position k:
    # Left side: constants + all variables before it
    # Middle: [variable_k, rhs]
    # Right side: all variables after it

    for i, var in enumerate(variables):
        left_vars = variables[:i]
        right_vars = variables[i+1:]

        # Build the new term:
        # Start with all constants:
        new_term_parts = constants[:]
        # Add left variables
        new_term_parts += left_vars
        # Insert the commutator for this variable
        new_term_parts.append(f"[{var}, {rhs}]")
        # Add right variables
        new_term_parts += right_vars

        new_term = '*'.join(new_term_parts)
        terms.append(new_term)

    # Combine terms with ' + '
    result = ' + '.join(terms)

    return result


leibniz_all_comms = [apply_leibniz_rule(all_comms_to_compute[i]) for i in range(len(all_comms_to_compute))]
leib_z_all_comms = [apply_leibniz_rule(all_z_comms_to_compute[i]) for i in range(len(all_z_comms_to_compute))]


def top_level_split_sign(expr):
    """
    Splits the input expression at top-level '+' or '-' signs,
    preserving the sign for each resulting piece. (A plus or minus
    is only considered a delimiter if it occurs at level 0, i.e.
    outside any parentheses.)

    For example:
      Input: "(n2+1)*Y3*Z3-n3*Y2*Z2"
      Output: [('+', "(n2+1)*Y3*Z3"), ('-', "n3*Y2*Z2")]
    """
    expr = expr.strip()
    result = []
    current = []
    level = 0
    i = 0

    # Set default sign based on first character:
    if expr and expr[0] in '+-':
        current_sign = expr[0]
        i = 1
    else:
        current_sign = '+'

    while i < len(expr):
        ch = expr[i]
        if ch == '(':
            level += 1
        elif ch == ')':
            level -= 1
        # If we are at level zero and see a plus or minus, that indicates a new term.
        if level == 0 and ch in '+-' and i != 0:
            term_str = ''.join(current).strip()
            if term_str:
                result.append((current_sign, term_str))
            current = []
            current_sign = ch
        else:
            current.append(ch)
        i += 1

    # Append the last term if any
    term_str = ''.join(current).strip()
    if term_str:
        result.append((current_sign, term_str))
    return [term[0] + term[1] for term in result]

def replace_minuses(equation):
    terms = top_level_split(equation)
    new_terms = []
    for term in terms:
        sub_terms_list = term.split('*')
        for var in sub_terms_list:
            if var[0] == '-':
                var_pos = sub_terms_list.index(var)
                sub_terms_list[var_pos] = var[1:]
                sub_terms_list.insert(var_pos, '(-1)')
        new_terms.append('*'.join(sub_terms_list))
    return ' + '.join(new_terms)

def delete_leading_1s(equation):
    terms = top_level_split(equation)
    new_terms = []
    for term in terms:
        if term[0] == '1':
            term = term[1:]
        elif term[1] == '1':
            term = term[0] + term[2:]
        new_terms.append(term)
    return ' + '.join(new_terms)

def replace_Z_commutators(expr):
    # Use top_level_split instead of split('+')
    terms = top_level_split(expr, '+')
    new_terms = []

    comm_pattern = re.compile(r'\[([^,]+),\s*a1[_^]([a-zA-Z])\*P1[_^]([a-zA-Z])\]')
    z_pattern = re.compile(r'^(Z)(\d+)(?:\^(.+))?$')

    for t in terms:
        term = t.strip()
        if not term:
            continue

        match = comm_pattern.search(term)
        if match:
            inside = match.group(1).strip()
            z_match = z_pattern.match(inside)
            if z_match:
                base = z_match.group(1)  # 'Z'
                number = z_match.group(2) # '1', '2', '3'
                exponent = z_match.group(3)
                if exponent is None:
                    exponent = "1"

                if number == '1':
                    # Z1 commutator -> remove entire term
                    continue
                elif number == '2':
                    # [Z2^n, a1_m*P1^m] -> n*Z2^(n-1)*Y3
                    new_exp = decrement_exponent(exponent)
                    replacement = f"{exponent}*Z2^{new_exp}*Y3"
                    new_term = term[:match.start()] + replacement + term[match.end():]
                    new_terms.append(new_term.strip())
                elif number == '3':
                    # [Z3^n, a1_m*P1^m] -> n*Z3^(n-1)*U2_a*P1^a
                    new_exp = decrement_exponent(exponent)
                    replacement = f"{exponent}*Z3^{new_exp}*U2_{find_unused_letter(t)}*P1^{find_unused_letter(t)}"
                    new_term = term[:match.start()] + replacement + term[match.end():]
                    new_terms.append(new_term.strip())
                else:
                    new_terms.append(term.strip())
            else:
                # Not a Z variable
                new_terms.append(term.strip())
        else:
            # no commutator
            new_terms.append(term.strip())

    result = ' + '.join(new_terms)
    return result

def apply_leibniz_to_y(expression):
    """
    Applies Leibniz' rule to commutators of the form [Y<id>^n, a1_m*P1^m].
    If n>1, expands into a sum of terms.
    If n=1, just leave it as [Y<id>, a1_m*P1^m].
    """

    terms = top_level_split(expression, '+')
    comm_pattern = re.compile(r'\[([^,]+),\s*a1[_^]([a-zA-Z])\*P1[_^]([a-zA-Z])\]')
    y_pattern = re.compile(r'^(Y\d+)(?:\^(\d+))?$')

    new_terms = []
    for t in terms:
        term = t.strip()
        if not term:
            continue

        match = comm_pattern.search(term)
        if match:
            inside = match.group(1).strip()
            y_match = y_pattern.match(inside)
            if y_match:
                y_base = y_match.group(1)  # e.g. 'Y1'
                exponent_str = y_match.group(2)
                if exponent_str is None:
                    exponent = 1
                else:
                    exponent = int(exponent_str)

                if exponent <= 1:
                    # Just [Y1, a1_m*P1^m], no expansion
                    new_terms.append(term)
                else:
                    # Apply Leibniz rule:
                    # [Y^n, A] = Σ_k=1^n Y^{n-k} [Y,A] Y^{k-1}
                    start, end = match.span()
                    left_part = term[:start]
                    right_part = term[end:]

                    expanded_terms = []
                    for k in range(1, exponent+1):
                        left_factors = []
                        if exponent - k > 0:
                            left_factors.append(f"{y_base}^{exponent - k}")

                        new_letter = find_unused_letter(t)
                        comm_str = f"[{y_base}, a1_{new_letter}*P1^{new_letter}]"

                        right_factors = []
                        if k - 1 > 0:
                            right_factors.append(f"{y_base}^{k - 1}")

                        # Construct each expanded term
                        # We'll rebuild carefully:
                        # Original term: left_part [Y_base^n, ...] right_part
                        # New term: left_part (left_factors + comm_str + right_factors) right_part
                        
                        # Clean up left_part and right_part
                        l = left_part.rstrip('*').strip()
                        r = right_part.lstrip('*').strip()

                        middle_factors = left_factors + [comm_str] + right_factors
                        middle_factors = [f for f in middle_factors if f]

                        # Combine all pieces with '*'
                        # We'll form a list and join by '*'
                        combined_factors = []
                        if l:
                            combined_factors.append(l.strip('*'))
                        if middle_factors:
                            combined_factors.append('*'.join(middle_factors))
                        if r:
                            combined_factors.append(r.strip('*'))

                        # Join all with '*', removing extra '*'
                        final_term = '*'.join(cf.strip('*') for cf in combined_factors if cf)
                        final_term = re.sub(r'\*+', '*', final_term)
                        final_term = final_term.strip('*')
                        expanded_terms.append(final_term)

                    new_terms.append(' + '.join(expanded_terms))
            else:
                # Not a Y commutator we care about
                new_terms.append(term)
        else:
            # no commutator in this term
            new_terms.append(term)

    # Post-processing: remove ^1 from Y variables
    # Replace Y<id>^1 with Y<id>
    result = ' + '.join(new_terms)
    result = re.sub(r'(Y\d+)\^1\b', r'\1', result)

    return result


all_leibniz_all_comms = [apply_leibniz_to_y(replace_Z_commutators(leibniz_all_comms[i])) for i in range(len(leibniz_all_comms))]
all_leib_z_all_comms = [apply_leibniz_to_y(replace_Z_commutators(leib_z_all_comms[i])) for i in range(len(leib_z_all_comms))]


def write_Y_powers_out(equation):
    terms = top_level_split(equation)
    resulting_terms = []
    for term in terms:
        sub_terms_list = term.split('*')
        for var in sub_terms_list:
            if var.startswith('Y') and '^' in var:
                Y_var, Y_power = var.split('^')
                new_Y_string = []
                for i in range(int(Y_power)):
                    new_Y_string.append(Y_var)
                new_Y_string = '*'.join(new_Y_string)
                sub_terms_list[sub_terms_list.index(var)] = new_Y_string
        resulting_terms.append('*'.join(sub_terms_list))
    return ' + '.join(resulting_terms)

def impose_y_commutation_rules(expr):
    """
    Generalizes the patterns:
      [Y1, a1_x*P1^x], [Y2, a1_x*P1^x], [Y3, a1_x*P1^x]
    to allow any single-letter index x with either
    subscript or superscript on a1 and P1, i.e.:
      a1_x*P1^x, a1^x*P1_x, a1^x*P1^x, a1_x*P1_x

    Then applies:
      - Y2 commutator terms are removed entirely.
      - Y1 commutator terms are replaced with 'P1_d*P2^d'.
      - Y3 commutator terms are left unchanged.
    """

    # Pattern that allows either ^ or _ on both a1 and P1, capturing the same single-letter index.
    # Explanation:
    #   \[Y1,\s*        matches literal `[Y1,` plus optional whitespace
    #   a1[_^]([a-zA-Z]) matches `a1` then either `_` or `^` followed by a single letter -> group(1)
    #   \*P1[_^]\1       matches `*P1`, then either `_` or `^` plus the same letter captured by group(1)
    #   \]               closing bracket
    y1_comm = r'\[Y1,\s*a1[_^]([a-zA-Z])\*P1[_^]\1\]'
    y2_comm = r'\[Y2,\s*a1[_^]([a-zA-Z])\*P1[_^]\1\]'
    y3_comm = r'\[Y3,\s*a1[_^]([a-zA-Z])\*P1[_^]\1\]'

    # Split the expression into terms at top-level '+'
    terms = top_level_split(expr, '+')
    new_terms = []

    for term in terms:
        # If the term matches [Y2, ...], remove it (skip).
        if re.search(y2_comm, term):
            continue

        # If the term matches [Y1, ...], replace with 'P1_d*P2^d'.
        new_letter = find_unused_letter(term)
        term = re.sub(y1_comm, f'P1_{new_letter}*P2^{new_letter}', term)

        # If the term matches [Y3, ...], do nothing (leave it).
        # So we neither skip nor replace Y3 commutators.

        new_terms.append(term)

    # If no terms survive, return '0'.
    if not new_terms:
        return '0'

    # Recombine with ' + '.
    result = ' + '.join(new_terms)

    # Clean up multiple '*' if they appear, and strip extra '*' if needed.
    result = re.sub(r'\*+', '*', result).strip('*')
    return result

def move_constant_factors_front(expr):
    terms = top_level_split(expr, '+')
    new_terms = []

    # Conditions for a factor to be considered a "constant":
    # 1. Parenthesized factor: ( ... ) with no '^'
    # 2. One of 'n1', 'n2', 'n3' exactly
    for term in terms:
        factors = [f.strip() for f in term.split('*') if f.strip()]

        constants = []
        non_constants = []
        for f in factors:
            if ((f.startswith('(') and f.endswith(')') and '^' not in f) or
                f in ['n1', 'n2', 'n3']):
                # This is a constant factor
                constants.append(f)
            else:
                non_constants.append(f)

        # Move all constants to the front
        reordered_factors = constants + non_constants
        new_term = '*'.join(reordered_factors)
        # Clean up any accidental multiple '*' 
        new_term = re.sub(r'\*+', '*', new_term).strip('*')
        new_terms.append(new_term)

    result = ' + '.join(new_terms)
    return result

def replace_Y3_with_U3P1(expr):
    """
    Replaces each instance of Y3 with U3_<x>*P1^<x>, using a UNIQUE 'x'
    index for each occurrence. 
    """

    # Keep track of letters used in the original expr so we don't reuse them.
    used_letters = set(re.findall(r'[a-zA-Z]', expr))

    def _replacement(_match):
        # For each Y3 matched, find an unused letter and insert it.
        for letter_to_try in "bcefghijklmnopqrstuvwxyz":
            if letter_to_try not in used_letters:
                used_letters.add(letter_to_try)
                return f'U3_{letter_to_try}*P1^{letter_to_try}'
        raise ValueError("Ran out of letters for substitution.")

    # Use the callback to handle each occurrence individually
    return re.sub(r'\bY3\b', _replacement, expr)

def replace_Y1_with_U1P2(expr):
    """
    Replaces each instance of Y1 with U1_<x>*P2^<x>, using a UNIQUE 'x'
    index for each occurrence. 
    """

    # Keep track of letters used in the original expr so we don't reuse them.
    used_letters = set(re.findall(r'[a-zA-Z]', expr))

    def _replacement(_match):
        # For each Y3 matched, find an unused letter and insert it.
        for letter_to_try in "bcefghijklmnopqrstuvwxyz":
            if letter_to_try not in used_letters:
                used_letters.add(letter_to_try)
                return f'U1_{letter_to_try}*P2^{letter_to_try}'
        raise ValueError("Ran out of letters for substitution.")

    # Use the callback to handle each occurrence individually
    return re.sub(r'\bY1\b', _replacement, expr)

def replace_Y2_with_U2P3(expr):
    """
    Replaces each instance of Y2 with U2_<x>*P3^<x>, using a UNIQUE 'x'
    index for each occurrence. 
    """

    # Keep track of letters used in the original expr so we don't reuse them.
    used_letters = set(re.findall(r'[a-zA-Z]', expr))

    def _replacement(_match):
        # For each Y2 matched, find an unused letter and insert it.
        for letter_to_try in "bcefghijklmnopqrstuvwxyz":
            if letter_to_try not in used_letters:
                used_letters.add(letter_to_try)
                return f'U2_{letter_to_try}*P3^{letter_to_try}'
        raise ValueError("Ran out of letters for substitution.")

    # Use the callback to handle each occurrence individually
    return re.sub(r'\bY2\b', _replacement, expr)

def replace_Z1_with_U2U3(expr):
    """
    Replaces each instance of Z1 with U2_<x>*U3^<x>, using a UNIQUE 'x'
    index for each occurrence. 
    """

    # Keep track of letters used in the original expr so we don't reuse them.
    used_letters = set(re.findall(r'[a-zA-Z]', expr))

    def _replacement(_match):
        # For each Y2 matched, find an unused letter and insert it.
        for letter_to_try in "bcefghijklmnopqrstuvwxyz":
            if letter_to_try not in used_letters:
                used_letters.add(letter_to_try)
                return f'U2_{letter_to_try}*U3^{letter_to_try}'
        raise ValueError("Ran out of letters for substitution.")

    # Use the callback to handle each occurrence individually
    return re.sub(r'\bZ1\b', _replacement, expr)

def replace_Z2_with_U1U3(expr):
    """
    Replaces each instance of Z2 with U1_<x>*U3^<x>, using a UNIQUE 'x'
    index for each occurrence. 
    """

    # Keep track of letters used in the original expr so we don't reuse them.
    used_letters = set(re.findall(r'[a-zA-Z]', expr))

    def _replacement(_match):
        # For each Z2 matched, find an unused letter and insert it.
        for letter_to_try in "bcefghijklmnopqrstuvwxyz":
            if letter_to_try not in used_letters:
                used_letters.add(letter_to_try)
                return f'U1_{letter_to_try}*U3^{letter_to_try}'
        raise ValueError("Ran out of letters for substitution.")

    # Use the callback to handle each occurrence individually
    return re.sub(r'\bZ2\b', _replacement, expr)

def replace_Z3_with_U1U2(expr):
    """
    Replaces each instance of Z3 with U1_<x>*U2^<x>, using a UNIQUE 'x'
    index for each occurrence. 
    """

    # Keep track of letters used in the original expr so we don't reuse them.
    used_letters = set(re.findall(r'[a-zA-Z]', expr))

    def _replacement(_match):
        # For each Z2 matched, find an unused letter and insert it.
        for letter_to_try in "bcefghijklmnopqrstuvwxyz":
            if letter_to_try not in used_letters:
                used_letters.add(letter_to_try)
                return f'U1_{letter_to_try}*U2^{letter_to_try}'
        raise ValueError("Ran out of letters for substitution.")

    # Use the callback to handle each occurrence individually
    return re.sub(r'\bZ3\b', _replacement, expr)

def replace_all_Ys(expr):
    return replace_Y3_with_U3P1(replace_Y2_with_U2P3(replace_Y1_with_U1P2(expr)))

def replace_all_YZs(expr):
    return replace_Y1_with_U1P2(replace_Y2_with_U2P3(replace_Y3_with_U3P1(replace_Z1_with_U2U3(replace_Z2_with_U1U3(replace_Z3_with_U1U2(expr))))))

def replace_Y3_commutator(expr):
    """
    For each occurrence of:
        [Y3, a1_x*P1^y]   or   [Y3, a1^x*P1_y]   etc.
    we replace it with:
        U3_<new_letter>*a1_x*[P1^<new_letter>, P1^y]
    ensuring that each occurrence uses a *unique* <new_letter>.

    Subscript/superscript is allowed in either place for a1 and P1.
    The letters x,y come from the original expression, and <new_letter>
    is a fresh letter chosen from 'abcdefghijklmnopqrstuvwxyz'
    that is not yet in the expression (or used by previous matches).
    """

    # Regex:
    # [Y3, a1_ or a1^ (group(1)), then *P1_ or P1^ (group(2))]
    pattern = r'\[Y3,\s*a1[_^]([a-zA-Z])\*P1[_^]([a-zA-Z])\]'

    # Collect all letters that appear in the expression
    # so we don't reuse them for the new index.
    used_letters = set(re.findall(r'[a-zA-Z]', expr))

    # Callback function that handles each match:
    def _replacement(match):
        # Original letters from the expression
        x = match.group(1)
        y = match.group(2)

        # Find the first unused letter
        for candidate in "abcdefghijklmnopqrstuvwxyz":
            if candidate not in used_letters:
                used_letters.add(candidate)
                # Build the replacement using x, y, and the newly chosen letter
                return f'U3_{candidate}*a1_{x}*[P1^{candidate},P1^{y}]'

        raise ValueError("Ran out of letters for substitution.")

    # Use the callback for each commutator match
    new_expr = re.sub(pattern, _replacement, expr)
    return new_expr

def replace_P1_commutator(expr):
    terms = top_level_split(expr)
    expanded_terms = []
    for term in terms:
        if '[P1' in term:
            # break into individual variables by * signs
            # this will be the main list which will be rearranged as we go
            sub_terms_list = term.split('*')
            # want to move all the U1, Z2, Z3 left of the leftmost P1 before we can commute the P1s
            # get all the U3, Z2, Z1 variables in a list - but only if they are on the right of the [P1,P1]
            # find [P1,P1] location
            for var in sub_terms_list:
                if var.startswith('[P1'):
                    init_P1_comm_pos = sub_terms_list.index(var)
            # then construct list of terms on the right of this
            U1_Z2_Z3_list = [var for var in sub_terms_list if (var.startswith('U3') or var.startswith('U2') or var.startswith('U1') or var.startswith('Z2') or var.startswith('Z3')) and sub_terms_list.index(var) > init_P1_comm_pos]
            # get the list without these variables
            remaining_list = [var for var in sub_terms_list if var not in U1_Z2_Z3_list]
            # get the position of the leftmost P1 in the remaining list
            for var in remaining_list:
                if var.startswith('[P1'):
                    leftmost_P1 = var
                    leftmost_P1_pos = remaining_list.index(var)
                    # now sequentially insert the U1_Z3_Z3_list on the left of the leftmost P1
                    for U1_Z2_Z3 in U1_Z2_Z3_list:
                        remaining_list.insert(leftmost_P1_pos, U1_Z2_Z3)
                    # assign rearranged list to our main sub_terms_list
                    sub_terms_list = remaining_list
                    break
                else:
                    continue
                
            contains_P1_commutator = False
            for P1_comm in sub_terms_list:
                if P1_comm.startswith('[P1'):
                    contains_P1_commutator = True
                    left_P1_commutator_index = P1_comm[3:5]
                    right_P1_commutator_index = P1_comm[8:10]
                    P1_comm_pos = sub_terms_list.index(P1_comm)
                    P1s_on_right = [var for var in sub_terms_list if var.startswith('P1') and sub_terms_list.index(var) > P1_comm_pos]
                    for P1_r in P1s_on_right:
                        # create copy to manipulate terms without altering master list
                        sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*P1_a*g_cb
                        sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*P1_b*g_ca
                        # construct the first term using the first copy
                        # replace P1_r in the copy list with P1_{index from left P1 in commutator}
                        sub_terms_list_copy_1[sub_terms_list_copy_1.index(P1_r)] = f"P1{left_P1_commutator_index}"
                        # then delete the commutator from the list
                        P1_comm_pos = sub_terms_list_copy_1.index(P1_comm)
                        sub_terms_list_copy_1.pop(P1_comm_pos)
                        # then in this term, we will also have a metric to contract: g_{P1_r index}_{index from right P1 in commutator}
                        # so find another term with one of these indices (choose P1_r cos why not) and change its letter index accordingly
                        for var in sub_terms_list_copy_1:
                            # check if index is same as P1_r -> this var will be contracted with the metric
                            if var[-1] == P1_r[-1]:
                                # find the index of this var in the list copy
                                sub_list_var_index = sub_terms_list_copy_1.index(var)
                                # if the P1_r index is being contracted with var, then the remaining index in the index of the
                                # P1 in the RHS of the commutator
                                sub_terms_list_copy_1[sub_list_var_index] = f"{var[:-2]}{right_P1_commutator_index}"
                                break
                        # now add in the 1/l^2 term, accounting for possible minus signs
                        if sub_terms_list_copy_1[0][0] != '-':
                            sub_terms_list_copy_1 = ['(1/l^2)'] + sub_terms_list_copy_1
                        else:
                            sub_terms_list_copy_1 = ['-(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                        expanded_terms.append('*'.join(sub_terms_list_copy_1))

                        # now repeat the process but for the second term; same structure, except for a sign and swap of indices
                        # construct the second term using the second copy
                        # replace the term in the LHS of the commutator in the copy list with P1_{index from right P3 in commutator}
                        sub_terms_list_copy_2[sub_terms_list_copy_2.index(P1_r)] = f"P1{right_P1_commutator_index}"
                        # delete the commutator
                        P1_comm_pos = sub_terms_list_copy_2.index(P1_comm)
                        sub_terms_list_copy_2.pop(P1_comm_pos)
                        # then in this term, we will also have a metric to contract: g_{P3_r index}_{index from left P3 in commutator}
                        # so find another term with one of these indices (choose P3_r cos why not) and change its letter index accordingly
                        for var in sub_terms_list_copy_2:
                            # check if index is same as P3_r -> this var will be contracted with the metric
                            if var[-1] == P1_r[-1]:
                                # find the index of this var in the list copy
                                sub_list_var_index = sub_terms_list_copy_2.index(var)
                                # if the P3_r index is being contracted with var, then the remaining index in the index of the
                                # P3 in the LHS of the commutator
                                sub_terms_list_copy_2[sub_list_var_index] = f"{var[:-2]}{left_P1_commutator_index}"
                                break
                        # now add in the 1/l^2 term, accounting for possible minus signs
                        if sub_terms_list_copy_2[0][0] != '-':
                            sub_terms_list_copy_2 = ['-(1/l^2)'] + sub_terms_list_copy_2
                        else:
                            sub_terms_list_copy_2 = ['(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                        expanded_terms.append('*'.join(sub_terms_list_copy_2))
                    # now we have to add in the riemann tensor operator that acts on the fields themselves.
                    # again, make two separate terms
                    sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*a1_b*U1_a
                    sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*a1_a*U1_b
                    # the indices are now only taken from P_var and P_var_on_left
                    # replace the commutator term with aU term
                    P1_commutator_index = sub_terms_list_copy_1.index(P1_comm)
                    sub_terms_list_copy_1[P1_commutator_index] = f"a1{left_P1_commutator_index}*U1{right_P1_commutator_index}"
                    # then add in the 1/l^2 term, accounting for signs
                    if sub_terms_list_copy_1[0][0] != '-':
                        sub_terms_list_copy_1 = ['-(1/l^2)'] + sub_terms_list_copy_1
                    else:
                        sub_terms_list_copy_1 = ['(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                    expanded_terms.append('*'.join(sub_terms_list_copy_1))

                    # then the second term
                    # replace the commutator term with aU term
                    P1_commutator_index = sub_terms_list_copy_2.index(P1_comm)
                    sub_terms_list_copy_2[P1_commutator_index] = f"a1{right_P1_commutator_index}*U1{left_P1_commutator_index}"
                    # then add in the 1/l^2 term, accounting for signs
                    if sub_terms_list_copy_2[0][0] != '-':
                        sub_terms_list_copy_2 = ['(1/l^2)'] + sub_terms_list_copy_2
                    else:
                        sub_terms_list_copy_2 = ['-(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
        else:  
            expanded_terms.append(term)
                
    
    return ' + '.join(expanded_terms)


rearrange_y_comms_all = [' + '.join([replace_P1_commutator(replace_Y2_with_U2P3(replace_Y1_with_U1P2(replace_Y3_with_U3P1(replace_Y3_commutator(move_constant_factors_front(impose_y_commutation_rules(write_Y_powers_out(term)))))))) for term in top_level_split(all_leibniz_all_comms[0])])]
rearrange_y_comms_all_z = [replace_P1_commutator(replace_Y2_with_U2P3(replace_Y1_with_U1P2(replace_Y3_with_U3P1(replace_Y3_commutator(move_constant_factors_front(impose_y_commutation_rules(write_Y_powers_out(all_leib_z_all_comms[i])))))))) for i in range(len(all_leib_z_all_comms))]


def move_a1_left(equation, drop_left_a1_terms=True):
    # break the equation down into individual terms by top level '+' symbols
    equation = replace_minuses(equation)
    terms = top_level_split(equation)
    expanded_terms = []
    for term in terms:
        # break into individual variables by * signs
        sub_terms_list = term.split('*')
        # pick out the list of relevant objects (in order) - the a1s, any explicit U1s, and Z2s and Z3s - alongside their positions in the original list
        a1s_Zs_Us_ordered = [(sub_terms_list[i], i) for i in range(len(sub_terms_list)) if sub_terms_list[i].startswith('U1') or sub_terms_list[i].startswith('Z3') or sub_terms_list[i].startswith('Z2') or sub_terms_list[i].startswith('a1')]
        # check for explicit a1s present in the filtered list - if not, we don't need to do anything
        if any(entry[0].startswith('a1') for entry in a1s_Zs_Us_ordered):
            # pick out just the terms, not the indices, for easier manipulation
            terms_to_combine = [a1s_Zs_Us_ordered[i][0] for i in range(len(a1s_Zs_Us_ordered))]
            # now look for the specific a1 in the list
            for var in terms_to_combine:
                if var.startswith('a1'):
                    # if the a1 is not on the left of all the operators i.e. its index is not zero, or all operators on the left are also a1s,
                    #  then we need to perform operations for rearrangement
                    all_operators_on_left = [var_1 for var_1 in terms_to_combine if terms_to_combine.index(var_1) < terms_to_combine.index(var)]
                    all_as_on_left = True
                    for operator in all_operators_on_left:
                        if not operator.startswith('a'):
                            all_as_on_left = False
                    while terms_to_combine.index(var) != 0 and not all_as_on_left:
                        # store the a1 index
                        a1_index = terms_to_combine.index(var)
                        # look at the term to the left of the a1 - need to commute it through this
                        var_on_left = terms_to_combine[a1_index-1]
                        # first check if term on the left is a Z3
                        if var_on_left.startswith('Z3'):
                            # swap order of terms in the reduced list
                            terms_to_combine[a1_index], terms_to_combine[a1_index-1] = terms_to_combine[a1_index-1], terms_to_combine[a1_index]
                            # update the original term list - swap a1 in its current position with the Z3 variable being commuted through
                            current_term_a1_index = sub_terms_list.index(var)
                            left_of_a1_term_index = sub_terms_list.index(var_on_left)
                            sub_terms_list[current_term_a1_index], sub_terms_list[left_of_a1_term_index] = sub_terms_list[left_of_a1_term_index], sub_terms_list[current_term_a1_index] 
                            # construct commuted term to add to expanded terms list
                            # first check if there is a power
                            # extract the index on the a1
                            a_lett_index = var[-2:]
                            if '^' in var_on_left:
                                # extract the power of Z3
                                Z_power = var_on_left.split('^')[-1]
                                # need to reduce the z power since [Z3^n, a1_b] = n*Z3^(n-1)*U2_b
                                reduced_Z_power = decrement_exponent(Z_power)
                                commutator_expr = f"Z3^{reduced_Z_power}*U2{a_lett_index}".strip()
                                # construct new term with commuted expression. Pull constant factor to the left. Create copy to avoid altering original
                                new_sub_terms_list = sub_terms_list.copy()
                                # place commutator expression where the a1 originally was
                                new_sub_terms_list[sub_terms_list.index(var)] = commutator_expr
                                # remove the 'var_on_left' since this is accounted for by the commutator
                                new_sub_terms_list = [x for i, x in enumerate(new_sub_terms_list) if i != new_sub_terms_list.index(var_on_left)]
                                # then add the constant factor on the left
                                new_sub_terms_list = [f"{Z_power}".strip()] + new_sub_terms_list
                                # then construct final combined string and add to expanded terms list
                                expanded_terms.append('*'.join(new_sub_terms_list))
                            else:
                                # no power of Z3, so just use [Z3, a1_b] = U2_b
                                commutator_expr = f"U2{a_lett_index}"
                                # construct new term with commuted expression. Create copy to avoid altering original
                                new_sub_terms_list = sub_terms_list.copy()
                                # place commutator expression where the a1 originally was
                                new_sub_terms_list[sub_terms_list.index(var)] = commutator_expr
                                # remove the 'var_on_left' since this is accounted for by the commutator
                                new_sub_terms_list = [x for i, x in enumerate(new_sub_terms_list) if i != new_sub_terms_list.index(var_on_left)]
                                expanded_terms.append('*'.join(new_sub_terms_list))

                        # next check if term on the left is a Z2
                        if var_on_left.startswith('Z2'):
                            # swap order of terms in the reduced list
                            terms_to_combine[a1_index], terms_to_combine[a1_index-1] = terms_to_combine[a1_index-1], terms_to_combine[a1_index]
                            # update the original term list - swap a1 in its current position with the Z2 variable being commuted through
                            current_term_a1_index = sub_terms_list.index(var)
                            left_of_a1_term_index = sub_terms_list.index(var_on_left)
                            sub_terms_list[current_term_a1_index], sub_terms_list[left_of_a1_term_index] = sub_terms_list[left_of_a1_term_index], sub_terms_list[current_term_a1_index] 
                            # construct commuted term to add to expanded terms list
                            # first check if there is a power
                            # extract the index on the a1
                            a_lett_index = var[-2:]
                            if '^' in var_on_left:
                                # extract the power of Z2
                                Z_power = var_on_left.split('^')[-1]
                                # need to reduce the z power since [Z2^n, a1_b] = n*Z2^(n-1)*U3_b
                                reduced_Z_power = decrement_exponent(Z_power)
                                commutator_expr = f"Z2^{reduced_Z_power}*U3{a_lett_index}".strip()
                                # construct new term with commuted expression. Pull constant factor to the left. Create copy to avoid altering original
                                new_sub_terms_list = sub_terms_list.copy()
                                # place commutator expression where the a1 originally was
                                new_sub_terms_list[sub_terms_list.index(var)] = commutator_expr
                                # remove the 'var_on_left' since this is accounted for by the commutator
                                new_sub_terms_list = [x for i, x in enumerate(new_sub_terms_list) if i != new_sub_terms_list.index(var_on_left)]
                                # then add the constant factor on the left
                                new_sub_terms_list = [f"{Z_power}".strip()] + new_sub_terms_list
                                # then construct final combined string and add to expanded terms list
                                expanded_terms.append('*'.join(new_sub_terms_list))
                            else:
                                # no power of Z2, so just use [Z2, a1_b] = U3_b
                                commutator_expr = f"U3{a_lett_index}"
                                # construct new term with commuted expression. Create copy to avoid altering original
                                new_sub_terms_list = sub_terms_list.copy()
                                # place commutator expression where the a1 originally was
                                new_sub_terms_list[sub_terms_list.index(var)] = commutator_expr
                                # remove the 'var_on_left' since this is accounted for by the commutator
                                new_sub_terms_list = [x for i, x in enumerate(new_sub_terms_list) if i != new_sub_terms_list.index(var_on_left)]
                                expanded_terms.append('*'.join(new_sub_terms_list))

                        # next check if term on the left is a U1
                        if var_on_left.startswith('U1'):
                            # extract the index on the a1
                            a_lett_index = var[-2:]
                            # swap order of terms in the reduced list
                            terms_to_combine[a1_index], terms_to_combine[a1_index-1] = terms_to_combine[a1_index-1], terms_to_combine[a1_index]
                            # construct commuted term to add to expanded terms list
                            # no powers of U1 will occur, so can just treat as [U1_a, a1_b] = delta(a, b)
                            # but if indices are the same, i.e. [U1_a, a1^a], then we can just multiply the expression by d (the dimension)
                            # so we need to find the contracted pair in the original list, and the commutator term will just be the original but with
                            # the contracted pair of U1 now holding a1s index. e.g. U1_h*P2^h*a1_b -> a1_b*U1_h*P2^h + P2_b
                            var_contracted_with_U1 = [var for var in sub_terms_list if var != var_on_left and var[-1] == var_on_left[-1]][0]
                            if var_contracted_with_U1 == var:
                                # then U1 is being contracted with the a1, so just delete U1 and a1 and insert a factor of d at the 
                                # start of sub_terms_list
                                new_sub_terms_list = sub_terms_list.copy()
                                U1_index = new_sub_terms_list.index(var_on_left)
                                new_sub_terms_list.pop(U1_index)
                                a1_index = new_sub_terms_list.index(var)
                                new_sub_terms_list.pop(a1_index)
                                new_sub_terms_list.insert(0, 'd')
                                
                            else:
                                # U1 not contracted with a1, so construct normal term
                                commutator_expr = f"{var_contracted_with_U1.split('^')[0].split('_')[0]}{a_lett_index}"
                                # update the original term list - swap a1 in its current position with the U1 variable being commuted through
                                current_term_a1_index = sub_terms_list.index(var)
                                left_of_a1_term_index = sub_terms_list.index(var_on_left)
                                sub_terms_list[current_term_a1_index], sub_terms_list[left_of_a1_term_index] = sub_terms_list[left_of_a1_term_index], sub_terms_list[current_term_a1_index] 
                                # construct new term with commuted expression. Create copy to avoid altering original
                                new_sub_terms_list = sub_terms_list.copy()
                                # place commutator expression where the a1 originally was
                                new_sub_terms_list[sub_terms_list.index(var_contracted_with_U1)] = commutator_expr
                                # remove the 'var_on_left' since this is accounted for by the commutator, and the original variable contracted with U1 to avoid double counting
                                new_sub_terms_list = [x for i, x in enumerate(new_sub_terms_list) if (i != new_sub_terms_list.index(var_on_left) and x != var)]
                            expanded_terms.append('*'.join(new_sub_terms_list))
                        # update all_as_on_left
                        all_operators_on_left = [var_1 for var_1 in terms_to_combine if terms_to_combine.index(var_1) < terms_to_combine.index(var)]
                        all_as_on_left = True
                        for operator in all_operators_on_left:
                            if not operator.startswith('a'):
                                all_as_on_left = False
                        
            # Then construct the original term but with a1 all the way on the left, if the condition is set
            if not drop_left_a1_terms:
                expanded_terms.append('*'.join(sub_terms_list))
                                
        else:
            # if nothing needs to be done, just return the original term
            expanded_terms.append('*'.join(sub_terms_list))
            continue

    return pull_minus_signs_to_front(' + '.join(expanded_terms))

def move_a2_left(equation, drop_left_a2_terms=True):
    # break the equation down into individual terms by top level '+' symbols
    equation = replace_minuses(equation)
    terms = top_level_split(equation)
    expanded_terms = []
    for term in terms:
        # break into individual variables by * signs
        sub_terms_list = term.split('*')
        # pick out the list of relevant objects (in order) - the a1s, any explicit U1s, and Z2s and Z3s - alongside their positions in the original list
        a2s_Zs_Us_ordered = [(sub_terms_list[i], i) for i in range(len(sub_terms_list)) if sub_terms_list[i].startswith('U2') or sub_terms_list[i].startswith('Z3') or sub_terms_list[i].startswith('Z1') or sub_terms_list[i].startswith('a2')]
        # check for explicit a1s present in the filtered list - if not, we don't need to do anything
        if any(entry[0].startswith('a2') for entry in a2s_Zs_Us_ordered):
            # pick out just the terms, not the indices, for easier manipulation
            terms_to_combine = [a2s_Zs_Us_ordered[i][0] for i in range(len(a2s_Zs_Us_ordered))]
            # now look for the specific a2 in the list
            for var in terms_to_combine:
                if var.startswith('a2'):
                    # if the a2 is not on the left of all the operators i.e. its index is not zero, or all operators on the left are also as,
                    #  then we need to perform operations for rearrangement
                    all_operators_on_left = [var_1 for var_1 in terms_to_combine if terms_to_combine.index(var_1) < terms_to_combine.index(var)]
                    all_as_on_left = True
                    for operator in all_operators_on_left:
                        if not operator.startswith('a'):
                            all_as_on_left = False
                    while terms_to_combine.index(var) != 0 and not all_as_on_left:
                        # store the a1 index
                        a2_index = terms_to_combine.index(var)
                        # look at the term to the left of the a1 - need to commute it through this
                        var_on_left = terms_to_combine[a2_index-1]
                        # first check if term on the left is a Z3
                        if var_on_left.startswith('Z3'):
                            # swap order of terms in the reduced list
                            terms_to_combine[a2_index], terms_to_combine[a2_index-1] = terms_to_combine[a2_index-1], terms_to_combine[a2_index]
                            # update the original term list - swap a2 in its current position with the Z3 variable being commuted through
                            current_term_a2_index = sub_terms_list.index(var)
                            left_of_a2_term_index = sub_terms_list.index(var_on_left)
                            sub_terms_list[current_term_a2_index], sub_terms_list[left_of_a2_term_index] = sub_terms_list[left_of_a2_term_index], sub_terms_list[current_term_a2_index] 
                            # construct commuted term to add to expanded terms list
                            # first check if there is a power
                            # extract the index on the a2
                            a_lett_index = var[-2:]
                            if '^' in var_on_left:
                                # extract the power of Z3
                                Z_power = var_on_left.split('^')[-1]
                                # need to reduce the z power since [Z3^n, a2_b] = n*Z3^(n-1)*U1_b
                                reduced_Z_power = decrement_exponent(Z_power)
                                commutator_expr = f"Z3^{reduced_Z_power}*U1{a_lett_index}".strip()
                                # construct new term with commuted expression. Pull constant factor to the left. Create copy to avoid altering original
                                new_sub_terms_list = sub_terms_list.copy()
                                # place commutator expression where the a2 originally was
                                new_sub_terms_list[sub_terms_list.index(var)] = commutator_expr
                                # remove the 'var_on_left' since this is accounted for by the commutator
                                new_sub_terms_list = [x for i, x in enumerate(new_sub_terms_list) if i != new_sub_terms_list.index(var_on_left)]
                                # then add the constant factor on the left
                                new_sub_terms_list = [f"{Z_power}".strip()] + new_sub_terms_list
                                # then construct final combined string and add to expanded terms list
                                expanded_terms.append('*'.join(new_sub_terms_list))
                            else:
                                # no power of Z3, so just use [Z3, a2_b] = U1_b
                                commutator_expr = f"U1{a_lett_index}"
                                # construct new term with commuted expression. Create copy to avoid altering original
                                new_sub_terms_list = sub_terms_list.copy()
                                # place commutator expression where the a2 originally was
                                new_sub_terms_list[sub_terms_list.index(var)] = commutator_expr
                                # remove the 'var_on_left' since this is accounted for by the commutator
                                new_sub_terms_list = [x for i, x in enumerate(new_sub_terms_list) if i != new_sub_terms_list.index(var_on_left)]
                                expanded_terms.append('*'.join(new_sub_terms_list))

                        # next check if term on the left is a Z1
                        if var_on_left.startswith('Z1'):
                            # swap order of terms in the reduced list
                            terms_to_combine[a2_index], terms_to_combine[a2_index-1] = terms_to_combine[a2_index-1], terms_to_combine[a2_index]
                            # update the original term list - swap a2 in its current position with the Z1 variable being commuted through
                            current_term_a2_index = sub_terms_list.index(var)
                            left_of_a2_term_index = sub_terms_list.index(var_on_left)
                            sub_terms_list[current_term_a2_index], sub_terms_list[left_of_a2_term_index] = sub_terms_list[left_of_a2_term_index], sub_terms_list[current_term_a2_index] 
                            # construct commuted term to add to expanded terms list
                            # first check if there is a power
                            # extract the index on the a2
                            a_lett_index = var[-2:]
                            if '^' in var_on_left:
                                # extract the power of Z1
                                Z_power = var_on_left.split('^')[-1]
                                # need to reduce the z power since [Z1^n, a2_b] = n*Z1^(n-1)*U3_b
                                reduced_Z_power = decrement_exponent(Z_power)
                                commutator_expr = f"Z1^{reduced_Z_power}*U3{a_lett_index}".strip()
                                # construct new term with commuted expression. Pull constant factor to the left. Create copy to avoid altering original
                                new_sub_terms_list = sub_terms_list.copy()
                                # place commutator expression where the a2 originally was
                                new_sub_terms_list[sub_terms_list.index(var)] = commutator_expr
                                # remove the 'var_on_left' since this is accounted for by the commutator
                                new_sub_terms_list = [x for i, x in enumerate(new_sub_terms_list) if i != new_sub_terms_list.index(var_on_left)]
                                # then add the constant factor on the left
                                new_sub_terms_list = [f"{Z_power}".strip()] + new_sub_terms_list
                                # then construct final combined string and add to expanded terms list
                                expanded_terms.append('*'.join(new_sub_terms_list))
                            else:
                                # no power of Z1, so just use [Z1, a2_b] = U3_b
                                commutator_expr = f"U3{a_lett_index}"
                                # construct new term with commuted expression. Create copy to avoid altering original
                                new_sub_terms_list = sub_terms_list.copy()
                                # place commutator expression where the a2 originally was
                                new_sub_terms_list[sub_terms_list.index(var)] = commutator_expr
                                # remove the 'var_on_left' since this is accounted for by the commutator
                                new_sub_terms_list = [x for i, x in enumerate(new_sub_terms_list) if i != new_sub_terms_list.index(var_on_left)]
                                expanded_terms.append('*'.join(new_sub_terms_list))

                        # next check if term on the left is a U2
                        if var_on_left.startswith('U2'):
                            # extract the index on the a2
                            a_lett_index = var[-2:]
                            # swap order of terms in the reduced list
                            terms_to_combine[a2_index], terms_to_combine[a2_index-1] = terms_to_combine[a2_index-1], terms_to_combine[a2_index]
                            # construct commuted term to add to expanded terms list
                            # no powers of U2 will occur, so can just treat as [U2_a, a2_b] = delta(a, b)
                            # so we need to find the contracted pair in the original list, and the commutator term will just be the original but with
                            # the contracted pair of U2 now holding a2s index. e.g. U2_h*P2^h*a2_b -> a2_b*U2_h*P2^h + P2_b
                            var_contracted_with_U2 = [var for var in sub_terms_list if var != var_on_left and var[-1] == var_on_left[-1]][0]
                            if var_contracted_with_U2 == var:
                                # then U1 is being contracted with the a1, so just delete U1 and a1 and insert a factor of d at the 
                                # start of sub_terms_list
                                new_sub_terms_list = sub_terms_list.copy()
                                U2_index = new_sub_terms_list.index(var_on_left)
                                new_sub_terms_list.pop(U2_index)
                                a2_index = new_sub_terms_list.index(var)
                                new_sub_terms_list.pop(a2_index)
                                new_sub_terms_list.insert(0, 'd')
                                
                            else:
                                # U2 not contracted with a2, so construct normal term
                                commutator_expr = f"{var_contracted_with_U2.split('^')[0].split('_')[0]}{a_lett_index}"
                                # update the original term list - swap a1 in its current position with the U1 variable being commuted through
                                current_term_a1_index = sub_terms_list.index(var)
                                left_of_a1_term_index = sub_terms_list.index(var_on_left)
                                sub_terms_list[current_term_a1_index], sub_terms_list[left_of_a1_term_index] = sub_terms_list[left_of_a1_term_index], sub_terms_list[current_term_a1_index] 
                                # construct new term with commuted expression. Create copy to avoid altering original
                                new_sub_terms_list = sub_terms_list.copy()
                                # place commutator expression where the a2 originally was
                                new_sub_terms_list[sub_terms_list.index(var_contracted_with_U2)] = commutator_expr
                                # remove the 'var_on_left' since this is accounted for by the commutator, and the original variable contracted with U2 to avoid double counting
                                new_sub_terms_list = [x for i, x in enumerate(new_sub_terms_list) if (i != new_sub_terms_list.index(var_on_left) and x != var)]
                            expanded_terms.append('*'.join(new_sub_terms_list))
                        # update all_as_on_left
                        all_operators_on_left = [var_1 for var_1 in terms_to_combine if terms_to_combine.index(var_1) < terms_to_combine.index(var)]
                        all_as_on_left = True
                        for operator in all_operators_on_left:
                            if not operator.startswith('a'):
                                all_as_on_left = False
                        
            # Then construct the original term but with a1 all the way on the left, if the condition is set
            if not drop_left_a2_terms:
                expanded_terms.append('*'.join(sub_terms_list))
                                
        else:
            # if nothing needs to be done, just return the original term
            expanded_terms.append('*'.join(sub_terms_list))
            continue

    return pull_minus_signs_to_front(' + '.join(expanded_terms))

def move_a3_left(equation, drop_left_a3_terms=True):
    # break the equation down into individual terms by top level '+' symbols
    equation = replace_minuses(equation)
    terms = top_level_split(equation)
    expanded_terms = []
    for term in terms:
        # break into individual variables by * signs
        sub_terms_list = term.split('*')
        # pick out the list of relevant objects (in order) - the a3s, any explicit U3s, and Z2s and Z1s - alongside their positions in the original list
        a3s_Zs_Us_ordered = [(sub_terms_list[i], i) for i in range(len(sub_terms_list)) if sub_terms_list[i].startswith('U3') or sub_terms_list[i].startswith('Z2') or sub_terms_list[i].startswith('Z1') or sub_terms_list[i].startswith('a3')]
        # check for explicit a1s present in the filtered list - if not, we don't need to do anything
        if any(entry[0].startswith('a3') for entry in a3s_Zs_Us_ordered):
            # pick out just the terms, not the indices, for easier manipulation
            terms_to_combine = [a3s_Zs_Us_ordered[i][0] for i in range(len(a3s_Zs_Us_ordered))]
            # now look for the specific a3 in the list
            for var in terms_to_combine:
                if var.startswith('a3'):
                    # if the a3 is not on the left of all the operators i.e. its index is not zero, or all operators on the left are also as,
                    #  then we need to perform operations for rearrangement
                    all_operators_on_left = [var_1 for var_1 in terms_to_combine if terms_to_combine.index(var_1) < terms_to_combine.index(var)]
                    all_as_on_left = True
                    for operator in all_operators_on_left:
                        if not operator.startswith('a'):
                            all_as_on_left = False
                    while terms_to_combine.index(var) != 0 and not all_as_on_left:
                        # store the a3 index
                        a3_index = terms_to_combine.index(var)
                        # look at the term to the left of the a3 - need to commute it through this
                        var_on_left = terms_to_combine[a3_index-1]
                        # first check if term on the left is a Z2
                        if var_on_left.startswith('Z2'):
                            # swap order of terms in the reduced list
                            terms_to_combine[a3_index], terms_to_combine[a3_index-1] = terms_to_combine[a3_index-1], terms_to_combine[a3_index]
                            # update the original term list - swap a3 in its current position with the Z2 variable being commuted through
                            current_term_a3_index = sub_terms_list.index(var)
                            left_of_a3_term_index = sub_terms_list.index(var_on_left)
                            sub_terms_list[current_term_a3_index], sub_terms_list[left_of_a3_term_index] = sub_terms_list[left_of_a3_term_index], sub_terms_list[current_term_a3_index] 
                            # construct commuted term to add to expanded terms list
                            # first check if there is a power
                            # extract the index on the a2
                            a_lett_index = var[-2:]
                            if '^' in var_on_left:
                                # extract the power of Z2
                                Z_power = var_on_left.split('^')[-1]
                                # need to reduce the z power since [Z2^n, a3_b] = n*Z2^(n-1)*U1_b
                                reduced_Z_power = decrement_exponent(Z_power)
                                commutator_expr = f"Z2^{reduced_Z_power}*U1{a_lett_index}".strip()
                                # construct new term with commuted expression. Pull constant factor to the left. Create copy to avoid altering original
                                new_sub_terms_list = sub_terms_list.copy()
                                # place commutator expression where the a3 originally was
                                new_sub_terms_list[sub_terms_list.index(var)] = commutator_expr
                                # remove the 'var_on_left' since this is accounted for by the commutator
                                new_sub_terms_list = [x for i, x in enumerate(new_sub_terms_list) if i != new_sub_terms_list.index(var_on_left)]
                                # then add the constant factor on the left
                                new_sub_terms_list = [f"{Z_power}".strip()] + new_sub_terms_list
                                # then construct final combined string and add to expanded terms list
                                expanded_terms.append('*'.join(new_sub_terms_list))
                            else:
                                # no power of Z2, so just use [Z2, a3_b] = U1_b
                                commutator_expr = f"U1{a_lett_index}"
                                # construct new term with commuted expression. Create copy to avoid altering original
                                new_sub_terms_list = sub_terms_list.copy()
                                # place commutator expression where the a3 originally was
                                new_sub_terms_list[sub_terms_list.index(var)] = commutator_expr
                                # remove the 'var_on_left' since this is accounted for by the commutator
                                new_sub_terms_list = [x for i, x in enumerate(new_sub_terms_list) if i != new_sub_terms_list.index(var_on_left)]
                                expanded_terms.append('*'.join(new_sub_terms_list))

                        # next check if term on the left is a Z1
                        if var_on_left.startswith('Z1'):
                            # swap order of terms in the reduced list
                            terms_to_combine[a3_index], terms_to_combine[a3_index-1] = terms_to_combine[a3_index-1], terms_to_combine[a3_index]
                            # update the original term list - swap a3 in its current position with the Z1 variable being commuted through
                            current_term_a3_index = sub_terms_list.index(var)
                            left_of_a3_term_index = sub_terms_list.index(var_on_left)
                            sub_terms_list[current_term_a3_index], sub_terms_list[left_of_a3_term_index] = sub_terms_list[left_of_a3_term_index], sub_terms_list[current_term_a3_index] 
                            # construct commuted term to add to expanded terms list
                            # first check if there is a power
                            # extract the index on the a3
                            a_lett_index = var[-2:]
                            if '^' in var_on_left:
                                # extract the power of Z1
                                Z_power = var_on_left.split('^')[-1]
                                # need to reduce the z power since [Z1^n, a3_b] = n*Z1^(n-1)*U2_b
                                reduced_Z_power = decrement_exponent(Z_power)
                                commutator_expr = f"Z1^{reduced_Z_power}*U2{a_lett_index}".strip()
                                # construct new term with commuted expression. Pull constant factor to the left. Create copy to avoid altering original
                                new_sub_terms_list = sub_terms_list.copy()
                                # place commutator expression where the a3 originally was
                                new_sub_terms_list[sub_terms_list.index(var)] = commutator_expr
                                # remove the 'var_on_left' since this is accounted for by the commutator
                                new_sub_terms_list = [x for i, x in enumerate(new_sub_terms_list) if i != new_sub_terms_list.index(var_on_left)]
                                # then add the constant factor on the left
                                new_sub_terms_list = [f"{Z_power}".strip()] + new_sub_terms_list
                                # then construct final combined string and add to expanded terms list
                                expanded_terms.append('*'.join(new_sub_terms_list))
                            else:
                                # no power of Z1, so just use [Z1, a3_b] = U2_b
                                commutator_expr = f"U2{a_lett_index}"
                                # construct new term with commuted expression. Create copy to avoid altering original
                                new_sub_terms_list = sub_terms_list.copy()
                                # place commutator expression where the a3 originally was
                                new_sub_terms_list[sub_terms_list.index(var)] = commutator_expr
                                # remove the 'var_on_left' since this is accounted for by the commutator
                                new_sub_terms_list = [x for i, x in enumerate(new_sub_terms_list) if i != new_sub_terms_list.index(var_on_left)]
                                expanded_terms.append('*'.join(new_sub_terms_list))

                        # next check if term on the left is a U3
                        if var_on_left.startswith('U3'):
                            # extract the index on the a2
                            a_lett_index = var[-2:]
                            # swap order of terms in the reduced list
                            terms_to_combine[a3_index], terms_to_combine[a3_index-1] = terms_to_combine[a3_index-1], terms_to_combine[a3_index]
                            # construct commuted term to add to expanded terms list
                            # no powers of U3 will occur, so can just treat as [U3_a, a3_b] = delta(a, b)
                            # so we need to find the contracted pair in the original list, and the commutator term will just be the original but with
                            # the contracted pair of U3 now holding a3s index. e.g. U3_h*P2^h*a3_b -> a3_b*U3_h*P2^h + P2_b
                            var_contracted_with_U3 = [var for var in sub_terms_list if var != var_on_left and var[-1] == var_on_left[-1]][0]
                            if var_contracted_with_U3 == var:
                                # then U3 is being contracted with the a3, so just delete U3 and a3 and insert a factor of d at the 
                                # start of sub_terms_list
                                new_sub_terms_list = sub_terms_list.copy()
                                U3_index = new_sub_terms_list.index(var_on_left)
                                new_sub_terms_list.pop(U3_index)
                                a3_index = new_sub_terms_list.index(var)
                                new_sub_terms_list.pop(a3_index)
                                new_sub_terms_list.insert(0, 'd')
                                
                            else:
                                # U3 not contracted with a3, so construct normal term
                                commutator_expr = f"{var_contracted_with_U3.split('^')[0].split('_')[0]}{a_lett_index}"
                                # update the original term list - swap a3 in its current position with the U3 variable being commuted through
                                current_term_a1_index = sub_terms_list.index(var)
                                left_of_a1_term_index = sub_terms_list.index(var_on_left)
                                sub_terms_list[current_term_a1_index], sub_terms_list[left_of_a1_term_index] = sub_terms_list[left_of_a1_term_index], sub_terms_list[current_term_a1_index] 
                                # construct new term with commuted expression. Create copy to avoid altering original
                                new_sub_terms_list = sub_terms_list.copy()
                                # place commutator expression where the index being contracted originally was
                                new_sub_terms_list[sub_terms_list.index(var_contracted_with_U3)] = commutator_expr
                                # remove the 'var_on_left' since this is accounted for by the commutator, and the a variable
                                new_sub_terms_list = [x for i, x in enumerate(new_sub_terms_list) if (i != new_sub_terms_list.index(var_on_left) and x != var)]
                            expanded_terms.append('*'.join(new_sub_terms_list))
                        # update all_as_on_left
                        all_operators_on_left = [var_1 for var_1 in terms_to_combine if terms_to_combine.index(var_1) < terms_to_combine.index(var)]
                        all_as_on_left = True
                        for operator in all_operators_on_left:
                            if not operator.startswith('a'):
                                all_as_on_left = False
                        
            # Then construct the original term but with a1 all the way on the left, if the condition is set
            if not drop_left_a3_terms:
                expanded_terms.append('*'.join(sub_terms_list))
                                
        else:
            # if nothing needs to be done, just return the original term
            expanded_terms.append('*'.join(sub_terms_list))
            continue
    return pull_minus_signs_to_front(' + '.join(expanded_terms))

def move_all_a_left(equation, drop_left_a_terms=True):
    # shorcut to apply all a term rearrangements
    a1s_left = move_a1_left(equation, drop_left_a1_terms=drop_left_a_terms)
    a2s_left = move_a2_left(a1s_left, drop_left_a2_terms=drop_left_a_terms)
    a3s_left = move_a3_left(a2s_left, drop_left_a3_terms=drop_left_a_terms)
    return a3s_left

def reorder_un_terms(equation):
    equation = replace_minuses(equation)
    def get_letter(term):
        if '^' in term:
            return term.split('^')[1]
        elif '_' in term:
            return term.split('_')[1]
        return ''
    terms = top_level_split(equation)
    U_terms = ['U1', 'U2', 'U3']
    rearranged_terms = []
    # iterate for each possible U terms
    for term in terms:
        result = term.split('*')
        for U in U_terms:
            # Find and move U terms
            i = 0
            while i < len(result):
                # find all entries containing particular U instance
                if result[i].startswith(U):
                    # find the index letter so we can find its contracted partner
                    current_letter = get_letter(result[i])
                
                    # Find its contracted partner
                    for j in range(len(result)):
                        if get_letter(result[j]) == current_letter and j != i:
                            u_term = result.pop(i)
                            # Insert U term to desired position
                            if j > i:
                                # if the matched term is after the term being moved, the pop above will lower the matched term's index by one
                                result.insert(j-1, u_term)
                                # moving a term to the right of the position currently being checked pushes the index of the next term we want to check
                                # down by one, so we don't need to increment i - subtract one to cancel the plus one at the end
                                # but note that if the U is only moved up by one position, we don't want to check it again, otherwise we get an infinite loop
                                if j > i + 1:
                                    i -= 1
                            else:
                                # otherwise, the matched term's index is unchanged
                                # moving a term to the left of the position currently being checked does not affect the index of the next term we want
                                # to check, so leave i alone 
                                result.insert(j, u_term)
                            break
                i += 1
        rearranged_terms.append('*'.join(result))
    
    return pull_minus_signs_to_front(' + '.join(rearranged_terms))

def pull_non_canon_U1P3s_left(equation):
    # pulls the leftmost U1P3 to the left to integrate by parts and put in canonical form
    # break the equation down into individual terms by top level '+' symbols
    equation = replace_minuses(equation)
    terms = top_level_split(equation)
    expanded_terms = []
    for term in terms:
        # break into individual variables by * signs
        # this will be the main list which will be rearranged as we go
        sub_terms_list = term.split('*')
        # check over sub_terms_list for U1 variables
        U1P3_present = False
        for U_var in sub_terms_list:
            if U1P3_present:
                break
            elif U_var.startswith('U1'):
                # store the letter index and position index of the U1
                U1_lett_index = U_var[-1] 
                U1_pos_index = sub_terms_list.index(U_var)
                for P_var in sub_terms_list:
                    if P_var.startswith('P3') and P_var[-1] == U1_lett_index:
                        U1P3_present = True
                        # store the P3 index (letter index should be same as U1 if relevant)
                        P3_pos_index = sub_terms_list.index(P_var)
                        # we have found a contracted U1P3, so then we can do the rearrangement
                        # now we want to move all the U3, Z2, Z1 left of the leftmost P3 before we can commute the P3s
                        # get all the U3, Z2, Z1 variables in a list
                        U3_Z2_Z1_list = [var for var in sub_terms_list if var.startswith('U1') or var.startswith('U2') or var.startswith('U3') or var.startswith('Z2') or var.startswith('Z1') or var.startswith('Z3')]
                        # get the list without these variables
                        remaining_list = [var for var in sub_terms_list if var not in U3_Z2_Z1_list]
                        # get the position of the leftmost P3 in the remaining list
                        for var in remaining_list:
                            if var.startswith('P3'):
                                leftmost_P3 = var
                                leftmost_P3_pos = remaining_list.index(var)
                                break
                        # now sequentially insert the U1_Z3_Z3_list on the left of the leftmost P3
                        for U3_Z2_Z1 in U3_Z2_Z1_list:
                            remaining_list.insert(leftmost_P3_pos, U3_Z2_Z1)
                        # assign rearranged list to our main sub_terms_list
                        sub_terms_list = remaining_list
                        # now we can start commuting P3s. We make a list of all the P3s
                        P3_to_commute = [var for var in sub_terms_list if var.startswith('P3')]
                        # want our P3 in the U1P3 contraction to be moved left - if not in position zero of P3_to_commute, carry on
                        # rearranging 
                        while P3_to_commute.index(P_var) != 0:
                            # find P3 index in the to_commute list
                            P3_index_commute_list = P3_to_commute.index(P_var)
                            # find P3 index in the sub_terms_list
                            P3_index_sub_list = sub_terms_list.index(P_var)
                            # store the P3 variable on the next left to be commuted through
                            P3_var_on_left_commute_list = P3_to_commute[P3_index_commute_list-1]
                            # find this variable's position in the sub_terms_list
                            P3_var_on_left_sub_list = sub_terms_list.index(P3_var_on_left_commute_list)
                            # extract relevant information about the P3s being commuted (indices)
                            # used for constructing the indices in the riemann tensor
                            left_P3_commutator_index = P3_var_on_left_commute_list[-2:]
                            right_P3_commutator_index = P_var[-2:]
                            # create terms arising from the riemann tensor contracting with P3s on the right
                            P3s_on_right = [var for var in P3_to_commute if P3_to_commute.index(var) > P3_index_commute_list]
                            for P3_r in P3s_on_right:
                                # create copy to manipulate terms without altering master list
                                sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*P1_a*g_cb
                                sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*P1_b*g_ca
                                # construct the first term using the first copy
                                # replace P3_r in the copy list with P3_{index from left P3 in commutator}
                                sub_terms_list_copy_1[sub_terms_list_copy_1.index(P3_r)] = f"P3{left_P3_commutator_index}"
                                # then delete the P3 that was in the LHS of the commutator the copy list since this is accounted for in the commutator
                                left_P3_commutator_pos = sub_terms_list_copy_1.index(P3_var_on_left_commute_list)
                                sub_terms_list_copy_1.pop(left_P3_commutator_pos)
                                # then also delete the P3 that was in the RHS of the commutator from the list
                                right_P3_commutator_pos = sub_terms_list_copy_1.index(P_var)
                                sub_terms_list_copy_1.pop(right_P3_commutator_pos)
                                # then in this term, we will also have a metric to contract: g_{P3_r index}_{index from right P3 in commutator}
                                # so find another term with one of these indices (choose P3_r cos why not) and change its letter index accordingly
                                if P_var[-1] == P3_r[-1]:
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['d*(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))
                                else:
                                    for var in sub_terms_list_copy_1:
                                        # check if index is same as P3_r -> this var will be contracted with the metric
                                        if var[-1] == P3_r[-1]:
                                            # find the index of this var in the list copy
                                            sub_list_var_index = sub_terms_list_copy_1.index(var)
                                            # if the P3_r index is being contracted with var, then the remaining index in the index of the
                                            # P3 in the RHS of the commutator
                                            sub_terms_list_copy_1[sub_list_var_index] = f"{var[:-2]}{right_P3_commutator_index}"
                                            break
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))

                                # now repeat the process but for the second term; same structure, except for a sign and swap of indices
                                # construct the second term using the second copy
                                # delete the P3_r in the copy list since this is accounted for in the commutator
                                left_P3_commutator_pos = sub_terms_list_copy_2.index(P3_var_on_left_commute_list)
                                sub_terms_list_copy_2.pop(left_P3_commutator_pos)
                                # then also delete the P3 that was in the RHS of the commutator from the list
                                right_P3_commutator_pos = sub_terms_list_copy_2.index(P_var)
                                sub_terms_list_copy_2.pop(right_P3_commutator_pos)
                                # replace P3_r in the copy list with P3_{index from right P3 in commutator}
                                sub_terms_list_copy_2[sub_terms_list_copy_2.index(P3_r)] = f"P3{right_P3_commutator_index}"
                                # then in this term, we will also have a metric to contract: g_{P3_r index}_{index from left P3 in commutator}
                                # so find another term with one of these indices (choose P3_r cos why not) and change its letter index accordingly
                                if P3_var_on_left_commute_list[-1] == P3_r[-1]:
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-d*(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                                else:
                                    for var in sub_terms_list_copy_2:
                                        # check if index is same as P3_r -> this var will be contracted with the metric
                                        if var[-1] == P3_r[-1]:
                                            # find the index of this var in the list copy
                                            sub_list_var_index = sub_terms_list_copy_2.index(var)
                                            # if the P3_r index is being contracted with var, then the remaining index in the index of the
                                            # P3 in the LHS of the commutator
                                            sub_terms_list_copy_2[sub_list_var_index] = f"{var[:-2]}{left_P3_commutator_index}"
                                            break
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                            
                            # now we have to add in the riemann tensor operator that acts on the fields themselves.
                            # again, make two separate terms
                            sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*a1_b*U1_a
                            sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*a1_a*U1_b
                            # the indices are now only taken from P_var and P_var_on_left
                            # replace the LHS commutator P3 term with aU term
                            LHS_commutator_P3_index = sub_terms_list_copy_1.index(P3_var_on_left_commute_list)
                            sub_terms_list_copy_1[LHS_commutator_P3_index] = f"a3{left_P3_commutator_index}*U3{right_P3_commutator_index}"
                            # then delete the RHS commutator P3 term from the list
                            RHS_commutator_P3_index = sub_terms_list_copy_1.index(P_var)
                            sub_terms_list_copy_1.pop(RHS_commutator_P3_index)
                            # then add in the 1/l^2 term, accounting for signs
                            if sub_terms_list_copy_1[0][0] != '-':
                                sub_terms_list_copy_1 = ['-(1/l^2)'] + sub_terms_list_copy_1
                            else:
                                sub_terms_list_copy_1 = ['(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_1))

                            # then the second term
                            # replace the LHS commutator P3 term with aU term
                            LHS_commutator_P3_index = sub_terms_list_copy_2.index(P3_var_on_left_commute_list)
                            sub_terms_list_copy_2[LHS_commutator_P3_index] = f"a3{right_P3_commutator_index}*U3{left_P3_commutator_index}"
                            # then delete the RHS commutator P3 term from the list
                            RHS_commutator_P3_index = sub_terms_list_copy_2.index(P_var)
                            sub_terms_list_copy_2.pop(RHS_commutator_P3_index)
                            # then add in the 1/l^2 term, accounting for signs
                            if sub_terms_list_copy_2[0][0] != '-':
                                sub_terms_list_copy_2 = ['(1/l^2)'] + sub_terms_list_copy_2
                            else:
                                sub_terms_list_copy_2 = ['-(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # lastly swap the orders of P3s in both lists for the next iteration
                            sub_terms_list[P3_var_on_left_sub_list], sub_terms_list[P3_index_sub_list] = sub_terms_list[P3_index_sub_list], sub_terms_list[P3_var_on_left_sub_list]
                            P3_to_commute[P3_index_commute_list-1], P3_to_commute[P3_index_commute_list] = P3_to_commute[P3_index_commute_list], P3_to_commute[P3_index_commute_list-1]
                        # only want to do this for the leftmost U1P3
                        break
        if U1P3_present:
            # now that sub_terms_list has the variables in the desired order, we can replace U1P3 with -U1P2 + -U1P1
            # create copies again to construct desired terms
            sub_terms_list_Y1 = sub_terms_list.copy()
            sub_terms_list_Div1 = sub_terms_list.copy()
            # can just use P_var and U_var - in first list, replace the P_var with P2 (same index), and in the second, replace with P1 (same index)
            # but need the replaced terms to be on the left of all derivatives due to the partial integration
            P_var_index = sub_terms_list.index(P_var)
            leftmost_P_index = P_var_index
            for var in sub_terms_list_Y1:
                # find left most P index
                if var.startswith('P'):
                    leftmost_P_index = sub_terms_list_Y1.index(var)
                    break
            # insert the new terms on the left of the leftmost P_index
            sub_terms_list_Y1.insert(leftmost_P_index, f'P2{P_var[-2:]}')
            sub_terms_list_Div1.insert(leftmost_P_index, f'P1{P_var[-2:]}')
            # find P_var index again, and delete this from the list
            P_var_index = sub_terms_list_Y1.index(P_var)
            sub_terms_list_Y1.pop(P_var_index)
            sub_terms_list_Div1.pop(P_var_index)
            # then we need to flip the overall sign due to integrating by parts
            if sub_terms_list_Y1[0][0] == '-':
                sub_terms_list_Y1 = [sub_terms_list_Y1[0][1:]] + sub_terms_list_Y1[1:]
            else:
                sub_terms_list_Y1 = [f'-{sub_terms_list_Y1[0]}'] + sub_terms_list_Y1[1:]
            if sub_terms_list_Div1[0][0] == '-':
                sub_terms_list_Div1 = [sub_terms_list_Div1[0][1:]] + sub_terms_list_Div1[1:]
            else:
                sub_terms_list_Div1 = [f'-{sub_terms_list_Div1[0]}'] + sub_terms_list_Div1[1:]
            # then we can construct the terms by joining with *
            expanded_terms.append('*'.join(sub_terms_list_Y1))
            expanded_terms.append('*'.join(sub_terms_list_Div1))
        else:
            # just reconstruct the original term
            expanded_terms.append('*'.join(sub_terms_list))

    return pull_minus_signs_to_front(' + '.join(expanded_terms))

def pull_non_canon_U2P1s_left(equation):
    # pulls the leftmost U2P1 to the left to integrate by parts and put in canonical form
    # break the equation down into individual terms by top level '+' symbols
    equation = replace_minuses(equation)
    terms = top_level_split(equation)
    expanded_terms = []
    for term in terms:
        # break into individual variables by * signs
        # this will be the main list which will be rearranged as we go
        sub_terms_list = term.split('*')
        # check over sub_terms_list for U2 variables
        U2P1_present = False
        for U_var in sub_terms_list:
            if U2P1_present:
                break
            elif U_var.startswith('U2'):
                # store the letter index and position index of the U2
                U2_lett_index = U_var[-1] 
                U2_pos_index = sub_terms_list.index(U_var)
                for P_var in sub_terms_list:
                    if P_var.startswith('P1') and P_var[-1] == U2_lett_index:
                        U2P1_present = True
                        # store the P1 index (letter index should be same as U2 if relevant)
                        P1_pos_index = sub_terms_list.index(P_var)
                        # we have found a contracted U2P1, so then we can do the rearrangement
                        # now we want to move all the U1, Z2, Z3 left of the leftmost P1 before we can commute the P1s
                        # get all the U1, Z2, Z3 variables in a list
                        U1_Z2_Z3_list = [var for var in sub_terms_list if var.startswith('U1') or var.startswith('U2') or var.startswith('U3') or var.startswith('Z2') or var.startswith('Z1') or var.startswith('Z3')]
                        # get the list without these variables
                        remaining_list = [var for var in sub_terms_list if var not in U1_Z2_Z3_list]
                        # get the position of the leftmost P1 in the remaining list
                        for var in remaining_list:
                            if var.startswith('P1'):
                                leftmost_P1 = var
                                leftmost_P1_pos = remaining_list.index(var)
                                break
                        # now sequentially insert the U1_Z3_Z3_list on the left of the leftmost P3
                        for U1_Z2_Z3 in U1_Z2_Z3_list:
                            remaining_list.insert(leftmost_P1_pos, U1_Z2_Z3)
                        # assign rearranged list to our main sub_terms_list
                        sub_terms_list = remaining_list
                        # now we can start commuting P3s. We make a list of all the P3s
                        P1_to_commute = [var for var in sub_terms_list if var.startswith('P1')]
                        # want our P1 in the U2P1 contraction to be moved left - if not in position zero of P1_to_commute, carry on
                        # rearranging 
                        while P1_to_commute.index(P_var) != 0:
                            # find P1 index in the to_commute list
                            P1_index_commute_list = P1_to_commute.index(P_var)
                            # find P1 index in the sub_terms_list
                            P1_index_sub_list = sub_terms_list.index(P_var)
                            # store the P1 variable on the next left to be commuted through
                            P1_var_on_left_commute_list = P1_to_commute[P1_index_commute_list-1]
                            # find this variable's position in the sub_terms_list
                            P1_var_on_left_sub_list = sub_terms_list.index(P1_var_on_left_commute_list)
                            # extract relevant information about the P1s being commuted (indices)
                            # used for constructing the indices in the riemann tensor
                            left_P1_commutator_index = P1_var_on_left_commute_list[-2:]
                            right_P1_commutator_index = P_var[-2:]
                            # create terms arising from the riemann tensor contracting with P1s on the right
                            P1s_on_right = [var for var in P1_to_commute if P1_to_commute.index(var) > P1_index_commute_list]
                            for P1_r in P1s_on_right:
                                # create copy to manipulate terms without altering master list
                                sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*P1_a*g_cb
                                sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*P1_b*g_ca
                                # construct the first term using the first copy
                                # replace P1_r in the copy list with P1_{index from left P3 in commutator}
                                sub_terms_list_copy_1[sub_terms_list_copy_1.index(P1_r)] = f"P1{left_P1_commutator_index}"
                                # then delete the P3 that was in the LHS of the commutator the copy list since this is accounted for in the commutator
                                left_P1_commutator_pos = sub_terms_list_copy_1.index(P1_var_on_left_commute_list)
                                sub_terms_list_copy_1.pop(left_P1_commutator_pos)
                                # then also delete the P1 that was in the RHS of the commutator from the list
                                right_P1_commutator_pos = sub_terms_list_copy_1.index(P_var)
                                sub_terms_list_copy_1.pop(right_P1_commutator_pos)
                                # then in this term, we will also have a metric to contract: g_{P1_r index}_{index from right P1 in commutator}
                                # so find another term with one of these indices (choose P1_r cos why not) and change its letter index accordingly
                                if P_var[-1] == P1_r[-1]:
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['d*(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))
                                else:
                                    for var in sub_terms_list_copy_1:
                                        # check if index is same as P1_r -> this var will be contracted with the metric
                                        if var[-1] == P1_r[-1]:
                                            # find the index of this var in the list copy
                                            sub_list_var_index = sub_terms_list_copy_1.index(var)
                                            # if the P1_r index is being contracted with var, then the remaining index in the index of the
                                            # P1 in the RHS of the commutator
                                            sub_terms_list_copy_1[sub_list_var_index] = f"{var[:-2]}{right_P1_commutator_index}"
                                            break
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))

                                # now repeat the process but for the second term; same structure, except for a sign and swap of indices
                                # construct the second term using the second copy
                                # delete the P1_r in the copy list since this is accounted for in the commutator
                                left_P1_commutator_pos = sub_terms_list_copy_2.index(P1_var_on_left_commute_list)
                                sub_terms_list_copy_2.pop(left_P1_commutator_pos)
                                # then also delete the P1 that was in the RHS of the commutator from the list
                                right_P1_commutator_pos = sub_terms_list_copy_2.index(P_var)
                                sub_terms_list_copy_2.pop(right_P1_commutator_pos)
                                # replace P1_r in the copy list with P1_{index from right P1 in commutator}
                                sub_terms_list_copy_2[sub_terms_list_copy_2.index(P1_r)] = f"P1{right_P1_commutator_index}"
                                # then in this term, we will also have a metric to contract: g_{P1_r index}_{index from left P1 in commutator}
                                # so find another term with one of these indices (choose P1_r cos why not) and change its letter index accordingly
                                if P1_var_on_left_commute_list[-1] == P1_r[-1]:
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-d*(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                                else:
                                    for var in sub_terms_list_copy_2:
                                        # check if index is same as P3_r -> this var will be contracted with the metric
                                        if var[-1] == P1_r[-1]:
                                            # find the index of this var in the list copy
                                            sub_list_var_index = sub_terms_list_copy_2.index(var)
                                            # if the P1_r index is being contracted with var, then the remaining index in the index of the
                                            # P1 in the LHS of the commutator
                                            sub_terms_list_copy_2[sub_list_var_index] = f"{var[:-2]}{left_P1_commutator_index}"
                                            break
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                            
                            # now we have to add in the riemann tensor operator that acts on the fields themselves.
                            # again, make two separate terms
                            sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*a1_b*U1_a
                            sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*a1_a*U1_b
                            # the indices are now only taken from P_var and P_var_on_left
                            # replace the LHS commutator P1 term with aU term
                            LHS_commutator_P1_index = sub_terms_list_copy_1.index(P1_var_on_left_commute_list)
                            sub_terms_list_copy_1[LHS_commutator_P1_index] = f"a1{left_P1_commutator_index}*U1{right_P1_commutator_index}"
                            # then delete the RHS commutator P1 term from the list
                            RHS_commutator_P1_index = sub_terms_list_copy_1.index(P_var)
                            sub_terms_list_copy_1.pop(RHS_commutator_P1_index)
                            # then add in the 1/l^2 term, accounting for signs
                            if sub_terms_list_copy_1[0][0] != '-':
                                sub_terms_list_copy_1 = ['-(1/l^2)'] + sub_terms_list_copy_1
                            else:
                                sub_terms_list_copy_1 = ['(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_1))

                            # then the second term
                            # replace the LHS commutator P1 term with aU term
                            LHS_commutator_P1_index = sub_terms_list_copy_2.index(P1_var_on_left_commute_list)
                            sub_terms_list_copy_2[LHS_commutator_P1_index] = f"a1{right_P1_commutator_index}*U1{left_P1_commutator_index}"
                            # then delete the RHS commutator P1 term from the list
                            RHS_commutator_P1_index = sub_terms_list_copy_2.index(P_var)
                            sub_terms_list_copy_2.pop(RHS_commutator_P1_index)
                            # then add in the 1/l^2 term, accounting for signs
                            if sub_terms_list_copy_2[0][0] != '-':
                                sub_terms_list_copy_2 = ['(1/l^2)'] + sub_terms_list_copy_2
                            else:
                                sub_terms_list_copy_2 = ['-(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # lastly swap the orders of P3s in both lists for the next iteration
                            sub_terms_list[P1_var_on_left_sub_list], sub_terms_list[P1_index_sub_list] = sub_terms_list[P1_index_sub_list], sub_terms_list[P1_var_on_left_sub_list]
                            P1_to_commute[P1_index_commute_list-1], P1_to_commute[P1_index_commute_list] = P1_to_commute[P1_index_commute_list], P1_to_commute[P1_index_commute_list-1]
                        # only want to do this for the leftmost U1P3
                        break
        if U2P1_present:
            # now that sub_terms_list has the variables in the desired order, we can replace U2P1 with -U2P3 + -U2P2
            # create copies again to construct desired terms
            sub_terms_list_Y1 = sub_terms_list.copy()
            sub_terms_list_Div1 = sub_terms_list.copy()
            # can just use P_var and U_var - in first list, replace the P_var with P3 (same index), and in the second, replace with P2 (same index)
            # but need the replaced terms to be on the left of all derivatives due to the partial integration
            P_var_index = sub_terms_list.index(P_var)
            leftmost_P_index = P_var_index
            for var in sub_terms_list_Y1:
                # find left most P index
                if var.startswith('P'):
                    leftmost_P_index = sub_terms_list_Y1.index(var)
                    break
            # insert the new terms on the left of the leftmost P_index
            sub_terms_list_Y1.insert(leftmost_P_index, f'P3{P_var[-2:]}')
            sub_terms_list_Div1.insert(leftmost_P_index, f'P2{P_var[-2:]}')
            # find P_var index again, and delete this from the list
            P_var_index = sub_terms_list_Y1.index(P_var)
            sub_terms_list_Y1.pop(P_var_index)
            sub_terms_list_Div1.pop(P_var_index)
            # then we need to flip the overall sign due to integrating by parts
            if sub_terms_list_Y1[0][0] == '-':
                sub_terms_list_Y1 = [sub_terms_list_Y1[0][1:]] + sub_terms_list_Y1[1:]
            else:
                sub_terms_list_Y1 = [f'-{sub_terms_list_Y1[0]}'] + sub_terms_list_Y1[1:]
            if sub_terms_list_Div1[0][0] == '-':
                sub_terms_list_Div1 = [sub_terms_list_Div1[0][1:]] + sub_terms_list_Div1[1:]
            else:
                sub_terms_list_Div1 = [f'-{sub_terms_list_Div1[0]}'] + sub_terms_list_Div1[1:]
            # then we can construct the terms by joining with *
            expanded_terms.append('*'.join(sub_terms_list_Y1))
            expanded_terms.append('*'.join(sub_terms_list_Div1))
        else:
            # just reconstruct the original term
            expanded_terms.append('*'.join(sub_terms_list))

    return pull_minus_signs_to_front(' + '.join(expanded_terms))

def pull_non_canon_U3P2s_left(equation):
    # pulls the leftmost U3P2 to the left to integrate by parts and put in canonical form
    # break the equation down into individual terms by top level '+' symbols
    equation = replace_minuses(equation)
    terms = top_level_split(equation)
    expanded_terms = []
    for term in terms:
        # break into individual variables by * signs
        # this will be the main list which will be rearranged as we go
        sub_terms_list = term.split('*')
        # check over sub_terms_list for U2 variables
        U3P2_present = False
        for U_var in sub_terms_list:
            if U3P2_present:
                break
            elif U_var.startswith('U3'):
                # store the letter index and position index of the U3
                U3_lett_index = U_var[-1] 
                U3_pos_index = sub_terms_list.index(U_var)
                for P_var in sub_terms_list:
                    if P_var.startswith('P2') and P_var[-1] == U3_lett_index:
                        U3P2_present = True
                        # store the P2 index (letter index should be same as U3 if relevant)
                        P2_pos_index = sub_terms_list.index(P_var)
                        # we have found a contracted U3P2, so then we can do the rearrangement
                        # now we want to move all the U3, Z2, Z1 left of the leftmost P2 before we can commute the P2s
                        # get all the U3, Z2, Z1 variables in a list
                        U2_Z3_Z1_list = [var for var in sub_terms_list if var.startswith('U1') or var.startswith('U2') or var.startswith('U3') or var.startswith('Z2') or var.startswith('Z1') or var.startswith('Z3')]
                        # get the list without these variables
                        remaining_list = [var for var in sub_terms_list if var not in U2_Z3_Z1_list]
                        # get the position of the leftmost P2 in the remaining list
                        for var in remaining_list:
                            if var.startswith('P2'):
                                leftmost_P2 = var
                                leftmost_P2_pos = remaining_list.index(var)
                                break
                        # now sequentially insert the U2_Z2_Z1_list on the left of the leftmost P2
                        for U2_Z3_Z1 in U2_Z3_Z1_list:
                            remaining_list.insert(leftmost_P2_pos, U2_Z3_Z1)
                        # assign rearranged list to our main sub_terms_list
                        sub_terms_list = remaining_list
                        # now we can start commuting P2s. We make a list of all the P2s
                        P2_to_commute = [var for var in sub_terms_list if var.startswith('P2')]
                        # want our P2 in the U3P2 contraction to be moved left - if not in position zero of P2_to_commute, carry on
                        # rearranging 
                        while P2_to_commute.index(P_var) != 0:
                            # find P2 index in the to_commute list
                            P2_index_commute_list = P2_to_commute.index(P_var)
                            # find P2 index in the sub_terms_list
                            P2_index_sub_list = sub_terms_list.index(P_var)
                            # store the P2 variable on the next left to be commuted through
                            P2_var_on_left_commute_list = P2_to_commute[P2_index_commute_list-1]
                            # find this variable's position in the sub_terms_list
                            P2_var_on_left_sub_list = sub_terms_list.index(P2_var_on_left_commute_list)
                            # extract relevant information about the P2s being commuted (indices)
                            # used for constructing the indices in the riemann tensor
                            left_P2_commutator_index = P2_var_on_left_commute_list[-2:]
                            right_P2_commutator_index = P_var[-2:]
                            # create terms arising from the riemann tensor contracting with P2s on the right
                            P2s_on_right = [var for var in P2_to_commute if P2_to_commute.index(var) > P2_index_commute_list]
                            for P2_r in P2s_on_right:
                                # create copy to manipulate terms without altering master list
                                sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*P2_a*g_cb
                                sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*P2_b*g_ca
                                # construct the first term using the first copy
                                # replace P2_r in the copy list with P2_{index from left P2 in commutator}
                                sub_terms_list_copy_1[sub_terms_list_copy_1.index(P2_r)] = f"P2{left_P2_commutator_index}"
                                # then delete the P2 that was in the LHS of the commutator the copy list since this is accounted for in the commutator
                                left_P2_commutator_pos = sub_terms_list_copy_1.index(P2_var_on_left_commute_list)
                                sub_terms_list_copy_1.pop(left_P2_commutator_pos)
                                # then also delete the P2 that was in the RHS of the commutator from the list
                                right_P2_commutator_pos = sub_terms_list_copy_1.index(P_var)
                                sub_terms_list_copy_1.pop(right_P2_commutator_pos)
                                # then in this term, we will also have a metric to contract: g_{P2_r index}_{index from right P2 in commutator}
                                # so find another term with one of these indices (choose P2_r cos why not) and change its letter index accordingly
                                if P_var[-1] == P2_r[-1]:
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['d*(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))
                                else:
                                    for var in sub_terms_list_copy_1:
                                        # check if index is same as P2_r -> this var will be contracted with the metric
                                        if var[-1] == P2_r[-1]:
                                            # find the index of this var in the list copy
                                            sub_list_var_index = sub_terms_list_copy_1.index(var)
                                            # if the P2_r index is being contracted with var, then the remaining index in the index of the
                                            # P2 in the RHS of the commutator
                                            sub_terms_list_copy_1[sub_list_var_index] = f"{var[:-2]}{right_P2_commutator_index}"
                                            break
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))

                                # now repeat the process but for the second term; same structure, except for a sign and swap of indices
                                # construct the second term using the second copy
                                # delete the P2_r in the copy list since this is accounted for in the commutator
                                left_P1_commutator_pos = sub_terms_list_copy_2.index(P2_var_on_left_commute_list)
                                sub_terms_list_copy_2.pop(left_P2_commutator_pos)
                                # then also delete the P1 that was in the RHS of the commutator from the list
                                right_P2_commutator_pos = sub_terms_list_copy_2.index(P_var)
                                sub_terms_list_copy_2.pop(right_P2_commutator_pos)
                                # replace P2_r in the copy list with P2_{index from right P2 in commutator}
                                sub_terms_list_copy_2[sub_terms_list_copy_2.index(P2_r)] = f"P2{right_P2_commutator_index}"
                                # then in this term, we will also have a metric to contract: g_{P2_r index}_{index from left P2 in commutator}
                                # so find another term with one of these indices (choose P2_r cos why not) and change its letter index accordingly
                                if P2_var_on_left_commute_list[-1] == P2_r[-1]:
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-d*(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                                else:
                                    for var in sub_terms_list_copy_2:
                                        # check if index is same as P2_r -> this var will be contracted with the metric
                                        if var[-1] == P2_r[-1]:
                                            # find the index of this var in the list copy
                                            sub_list_var_index = sub_terms_list_copy_2.index(var)
                                            # if the P2_r index is being contracted with var, then the remaining index in the index of the
                                            # P2 in the LHS of the commutator
                                            sub_terms_list_copy_2[sub_list_var_index] = f"{var[:-2]}{left_P2_commutator_index}"
                                            break
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                            
                            # now we have to add in the riemann tensor operator that acts on the fields themselves.
                            # again, make two separate terms
                            sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*a2_b*U2_a
                            sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*a2_a*U2_b
                            # the indices are now only taken from P_var and P_var_on_left
                            # replace the LHS commutator P2 term with aU term
                            LHS_commutator_P2_index = sub_terms_list_copy_1.index(P2_var_on_left_commute_list)
                            sub_terms_list_copy_1[LHS_commutator_P2_index] = f"a2{left_P2_commutator_index}*U2{right_P2_commutator_index}"
                            # then delete the RHS commutator P2 term from the list
                            RHS_commutator_P2_index = sub_terms_list_copy_1.index(P_var)
                            sub_terms_list_copy_1.pop(RHS_commutator_P2_index)
                            # then add in the 1/l^2 term, accounting for signs
                            if sub_terms_list_copy_1[0][0] != '-':
                                sub_terms_list_copy_1 = ['-(1/l^2)'] + sub_terms_list_copy_1
                            else:
                                sub_terms_list_copy_1 = ['(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_1))

                            # then the second term
                            # replace the LHS commutator P2 term with aU term
                            LHS_commutator_P2_index = sub_terms_list_copy_2.index(P2_var_on_left_commute_list)
                            sub_terms_list_copy_2[LHS_commutator_P2_index] = f"a2{right_P2_commutator_index}*U2{left_P2_commutator_index}"
                            # then delete the RHS commutator P2 term from the list
                            RHS_commutator_P2_index = sub_terms_list_copy_2.index(P_var)
                            sub_terms_list_copy_2.pop(RHS_commutator_P2_index)
                            # then add in the 1/l^2 term, accounting for signs
                            if sub_terms_list_copy_2[0][0] != '-':
                                sub_terms_list_copy_2 = ['(1/l^2)'] + sub_terms_list_copy_2
                            else:
                                sub_terms_list_copy_2 = ['-(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # lastly swap the orders of P3s in both lists for the next iteration
                            sub_terms_list[P2_var_on_left_sub_list], sub_terms_list[P2_index_sub_list] = sub_terms_list[P2_index_sub_list], sub_terms_list[P2_var_on_left_sub_list]
                            P2_to_commute[P2_index_commute_list-1], P2_to_commute[P2_index_commute_list] = P2_to_commute[P2_index_commute_list], P2_to_commute[P2_index_commute_list-1]
                        # only want to do this for the leftmost U3P2
                        break
        if U3P2_present:
            # now that sub_terms_list has the variables in the desired order, we can replace U3P2 with -U3P1 + -U3P3
            # create copies again to construct desired terms
            sub_terms_list_Y3 = sub_terms_list.copy()
            sub_terms_list_Div3 = sub_terms_list.copy()
            # can just use P_var and U_var - in first list, replace the P_var with P1 (same index), and in the second, replace with P3 (same index)
            P_var_index = sub_terms_list.index(P_var)
            leftmost_P_index = P_var_index
            for var in sub_terms_list_Y3:
                # find left most P index
                if var.startswith('P'):
                    leftmost_P_index = sub_terms_list_Y3.index(var)
                    break
            # insert the new terms on the left of the leftmost P_index
            sub_terms_list_Y3.insert(leftmost_P_index, f'P1{P_var[-2:]}')
            sub_terms_list_Div3.insert(leftmost_P_index, f'P3{P_var[-2:]}')
            # find P_var index again, and delete this from the list
            P_var_index = sub_terms_list_Y3.index(P_var)
            sub_terms_list_Y3.pop(P_var_index)
            sub_terms_list_Div3.pop(P_var_index)
            # then we need to flip the overall sign due to integrating by parts
            if sub_terms_list_Y3[0][0] == '-':
                sub_terms_list_Y3 = [sub_terms_list_Y3[0][1:]] + sub_terms_list_Y3[1:]
            else:
                sub_terms_list_Y3 = [f'-{sub_terms_list_Y3[0]}'] + sub_terms_list_Y3[1:]
            if sub_terms_list_Div3[0][0] == '-':
                sub_terms_list_Div3 = [sub_terms_list_Div3[0][1:]] + sub_terms_list_Div3[1:]
            else:
                sub_terms_list_Div3 = [f'-{sub_terms_list_Div3[0]}'] + sub_terms_list_Div3[1:]
            # then we can construct the terms by joining with *
            expanded_terms.append('*'.join(sub_terms_list_Y3))
            expanded_terms.append('*'.join(sub_terms_list_Div3))
        else:
            # just reconstruct the original term
            expanded_terms.append('*'.join(sub_terms_list))

    return pull_minus_signs_to_front(' + '.join(expanded_terms))

def pull_all_non_canon_UP_left(equation):
    return move_all_a_left(pull_non_canon_U1P3s_left(move_all_a_left(pull_non_canon_U2P1s_left(move_all_a_left(pull_non_canon_U3P2s_left(move_all_a_left(equation)))))))

def pull_Div1_right(equation):
    # pulls the leftmost U1P1 to the right to discard it
    # break the equation down into individual terms by top level '+' symbols
    equation = replace_minuses(equation)
    terms = top_level_split(equation)
    expanded_terms = []
    for term in terms:
        # break into individual variables by * signs
        # this will be the main list which will be rearranged as we go
        sub_terms_list = term.split('*')
        # check over sub_terms_list for U1 variables
        U1P1_present = False
        for U_var in reversed(sub_terms_list):
            if U1P1_present:
                break
            elif U_var.startswith('U1'):
                # store the letter index and position index of the U1
                U1_lett_index = U_var[-1] 
                U1_pos_index = sub_terms_list.index(U_var)
                U1P1_present = False
                for P_var in sub_terms_list:
                    if P_var.startswith('P1') and P_var[-1] == U1_lett_index:
                        U1P1_present = True
                        # store the P1 index (letter index should be same as U1 if relevant)
                        P1_pos_index = sub_terms_list.index(P_var)
                        # we have found a contracted U1P1, so then we can do the rearrangement
                        # now we want to move all the U1, Z2, Z3 left of the leftmost P1 before we can commute the P1s
                        # get all the U1, Z2, Z3 variables in a list
                        U1_Z2_Z3_list = [var for var in sub_terms_list if var.startswith('U1') or var.startswith('U2') or var.startswith('U3') or var.startswith('Z2') or var.startswith('Z1') or var.startswith('Z3')]
                        # get the list without these variables
                        remaining_list = [var for var in sub_terms_list if var not in U1_Z2_Z3_list]
                        # get the position of the leftmost P1 in the remaining list
                        for var in remaining_list:
                            if var.startswith('P1'):
                                leftmost_P1 = var
                                leftmost_P1_pos = remaining_list.index(var)
                                break
                        # now sequentially insert the U1_Z2_Z3_list on the left of the leftmost P3
                        for U1_Z2_Z3 in U1_Z2_Z3_list:
                            remaining_list.insert(leftmost_P1_pos, U1_Z2_Z3)
                        # assign rearranged list to our main sub_terms_list
                        sub_terms_list = remaining_list
                        # now we can start commuting P1s. We make a list of all the P1s
                        P1_to_commute = [var for var in sub_terms_list if var.startswith('P1')]
                        # want our P1 in the U1P1 contraction to be moved right - if not in final position of P1_to_commute, carry on
                        # rearranging 
                        while P1_to_commute.index(P_var) != (len(P1_to_commute)-1):
                            # find P1 index in the to_commute list
                            P1_index_commute_list = P1_to_commute.index(P_var)
                            # find P1 index in the sub_terms_list
                            P1_index_sub_list = sub_terms_list.index(P_var)
                            # store the P1 variable on the next right to be commuted through
                            P1_var_on_right_commute_list = P1_to_commute[P1_index_commute_list+1]
                            # find this variable's position in the sub_terms_list
                            P1_var_on_right_sub_list = sub_terms_list.index(P1_var_on_right_commute_list)
                            # extract relevant information about the P3s being commuted (indices)
                            # used for constructing the indices in the riemann tensor
                            right_P1_commutator_index = P1_var_on_right_commute_list[-2:]
                            left_P1_commutator_index = P_var[-2:]
                            # create terms arising from the riemann tensor contracting with P1s on the right
                            P1s_on_right = [var for var in P1_to_commute if P1_to_commute.index(var) > P1_index_commute_list+1]
                            for P1_r in P1s_on_right:
                                # create copy to manipulate terms without altering master list
                                sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*P1_a*g_cb
                                sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*P1_b*g_ca
                                # construct the first term using the first copy
                                # replace P1_r in the copy list with P1_{index from left P1 in commutator}
                                sub_terms_list_copy_1[sub_terms_list_copy_1.index(P1_r)] = f"P1{left_P1_commutator_index}"
                                # then delete the P1 that was in the LHS of the commutator the copy list since this is accounted for in the commutator
                                left_P1_commutator_pos = sub_terms_list_copy_1.index(P1_var_on_right_commute_list)
                                sub_terms_list_copy_1.pop(left_P1_commutator_pos)
                                # then also delete the P1 that was in the RHS of the commutator from the list
                                right_P1_commutator_pos = sub_terms_list_copy_1.index(P_var)
                                sub_terms_list_copy_1.pop(right_P1_commutator_pos)
                                if P1_var_on_right_commute_list[-1] == P1_r[-1]:
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['d*(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))
                                else:
                                    # then in this term, we will also have a metric to contract: g_{P1_r index}_{index from right P1 in commutator}
                                    # so find another term with one of these indices (choose P1_r cos why not) and change its letter index accordingly
                                    for var in sub_terms_list_copy_1:
                                        # check if index is same as P1_r -> this var will be contracted with the metric
                                        if var[-1] == P1_r[-1]:
                                            # find the index of this var in the list copy
                                            sub_list_var_index = sub_terms_list_copy_1.index(var)
                                            # if the P1_r index is being contracted with var, then the remaining index in the index of the
                                            # P1 in the RHS of the commutator
                                            sub_terms_list_copy_1[sub_list_var_index] = f"{var[:-2]}{right_P1_commutator_index}"
                                            break
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))

                                # now repeat the process but for the second term; same structure, except for a sign and swap of indices
                                # construct the second term using the second copy
                                # delete the P1_r in the copy list since this is accounted for in the commutator
                                left_P1_commutator_pos = sub_terms_list_copy_2.index(P1_var_on_right_commute_list)
                                sub_terms_list_copy_2.pop(left_P1_commutator_pos)
                                # then also delete the P1 that was in the RHS of the commutator from the list
                                right_P3_commutator_pos = sub_terms_list_copy_2.index(P_var)
                                sub_terms_list_copy_2.pop(right_P1_commutator_pos)
                                # replace P1_r in the copy list with P1_{index from right P1 in commutator}
                                sub_terms_list_copy_2[sub_terms_list_copy_2.index(P1_r)] = f"P1{right_P1_commutator_index}"
                                # then in this term, we will also have a metric to contract: g_{P1_r index}_{index from left P1 in commutator}
                                # so find another term with one of these indices (choose P1_r cos why not) and change its letter index accordingly
                                if P_var[-1] == P1_r[-1]:
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-d*(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                                else:
                                    for var in sub_terms_list_copy_2:
                                        # check if index is same as P1_r -> this var will be contracted with the metric
                                        if var[-1] == P1_r[-1]:
                                            # find the index of this var in the list copy
                                            sub_list_var_index = sub_terms_list_copy_2.index(var)
                                            # if the P3_r index is being contracted with var, then the remaining index in the index of the
                                            # P3 in the LHS of the commutator
                                            sub_terms_list_copy_2[sub_list_var_index] = f"{var[:-2]}{left_P1_commutator_index}"
                                            break
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                            
                            # now we have to add in the riemann tensor operator that acts on the fields themselves.
                            # again, make two separate terms
                            sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*a1_b*U1_a
                            sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*a1_a*U1_b
                            # the indices are now only taken from P_var and P_var_on_left
                            # replace the LHS commutator P1 term with aU term
                            LHS_commutator_P1_index = sub_terms_list_copy_1.index(P1_var_on_right_commute_list)
                            sub_terms_list_copy_1[LHS_commutator_P1_index] = f"a1{left_P1_commutator_index}*U1{right_P1_commutator_index}"
                            # then delete the RHS commutator P1 term from the list
                            RHS_commutator_P1_index = sub_terms_list_copy_1.index(P_var)
                            sub_terms_list_copy_1.pop(RHS_commutator_P1_index)
                            # then add in the 1/l^2 term, accounting for signs
                            if sub_terms_list_copy_1[0][0] != '-':
                                sub_terms_list_copy_1 = ['-(1/l^2)'] + sub_terms_list_copy_1
                            else:
                                sub_terms_list_copy_1 = ['(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_1))

                            # then the second term
                            # replace the LHS commutator P1 term with aU term
                            LHS_commutator_P1_index = sub_terms_list_copy_2.index(P1_var_on_right_commute_list)
                            sub_terms_list_copy_2[LHS_commutator_P1_index] = f"a1{right_P1_commutator_index}*U1{left_P1_commutator_index}"
                            # then delete the RHS commutator P1 term from the list
                            RHS_commutator_P1_index = sub_terms_list_copy_2.index(P_var)
                            sub_terms_list_copy_2.pop(RHS_commutator_P1_index)
                            # then add in the 1/l^2 term, accounting for signs
                            if sub_terms_list_copy_2[0][0] != '-':
                                sub_terms_list_copy_2 = ['(1/l^2)'] + sub_terms_list_copy_2
                            else:
                                sub_terms_list_copy_2 = ['-(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # lastly swap the orders of P1s in both lists for the next iteration
                            sub_terms_list[P1_var_on_right_sub_list], sub_terms_list[P1_index_sub_list] = sub_terms_list[P1_index_sub_list], sub_terms_list[P1_var_on_right_sub_list]
                            P1_to_commute[P1_index_commute_list+1], P1_to_commute[P1_index_commute_list] = P1_to_commute[P1_index_commute_list], P1_to_commute[P1_index_commute_list+1]    
                        # only want to do this for the leftmost U1P3
                        break
        if not U1P1_present:
            # just reconstruct the original term
            expanded_terms.append('*'.join(sub_terms_list))
    
    return pull_minus_signs_to_front(' + '.join(expanded_terms))

def pull_Div2_right(equation):
    # pulls the leftmost U2P2 to the right to discard it
    # break the equation down into individual terms by top level '+' symbols
    equation = replace_minuses(equation)
    terms = top_level_split(equation)
    expanded_terms = []
    for term in terms:
        # break into individual variables by * signs
        # this will be the main list which will be rearranged as we go
        sub_terms_list = term.split('*')
        # check over sub_terms_list for U1 variables
        U2P2_present = False
        for U_var in reversed(sub_terms_list):
            if U2P2_present:
                break
            elif U_var.startswith('U2'):
                # store the letter index and position index of the U2
                U2_lett_index = U_var[-1] 
                U2_pos_index = sub_terms_list.index(U_var)
                for P_var in sub_terms_list:
                    if P_var.startswith('P2') and P_var[-1] == U2_lett_index:
                        U2P2_present = True
                        # store the P2 index (letter index should be same as U2 if relevant)
                        P2_pos_index = sub_terms_list.index(P_var)
                        # we have found a contracted U2P2, so then we can do the rearrangement
                        # now we want to move all the U2, Z1, Z3 left of the leftmost P2 before we can commute the P1s
                        # get all the U2, Z1, Z3 variables in a list
                        U2_Z1_Z3_list = [var for var in sub_terms_list if var.startswith('U1') or var.startswith('U2') or var.startswith('U3') or var.startswith('Z2') or var.startswith('Z1') or var.startswith('Z3')]
                        # get the list without these variables
                        remaining_list = [var for var in sub_terms_list if var not in U2_Z1_Z3_list]
                        # get the position of the leftmost P2 in the remaining list
                        for var in remaining_list:
                            if var.startswith('P2'):
                                leftmost_P2 = var
                                leftmost_P2_pos = remaining_list.index(var)
                                break
                        # now sequentially insert the U2_Z1_Z3_list on the left of the leftmost P2
                        for U2_Z1_Z3 in U2_Z1_Z3_list:
                            remaining_list.insert(leftmost_P2_pos, U2_Z1_Z3)
                        # assign rearranged list to our main sub_terms_list
                        sub_terms_list = remaining_list
                        # now we can start commuting P2s. We make a list of all the P2s
                        P2_to_commute = [var for var in sub_terms_list if var.startswith('P2')]
                        # want our P2 in the U2P2 contraction to be moved right - if not in final position of P2_to_commute, carry on
                        # rearranging 
                        while P2_to_commute.index(P_var) != (len(P2_to_commute)-1):
                            # find P2 index in the to_commute list
                            P2_index_commute_list = P2_to_commute.index(P_var)
                            # find P2 index in the sub_terms_list
                            P2_index_sub_list = sub_terms_list.index(P_var)
                            # store the P2 variable on the next right to be commuted through
                            P2_var_on_right_commute_list = P2_to_commute[P2_index_commute_list+1]
                            # find this variable's position in the sub_terms_list
                            P2_var_on_right_sub_list = sub_terms_list.index(P2_var_on_right_commute_list)
                            # extract relevant information about the P2s being commuted (indices)
                            # used for constructing the indices in the riemann tensor
                            right_P2_commutator_index = P2_var_on_right_commute_list[-2:]
                            left_P2_commutator_index = P_var[-2:]
                            # create terms arising from the riemann tensor contracting with P2s on the right
                            P2s_on_right = [var for var in P2_to_commute if P2_to_commute.index(var) > P2_index_commute_list+1]
                            for P2_r in P2s_on_right:
                                # create copy to manipulate terms without altering master list
                                sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*P2_a*g_cb
                                sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*P2_b*g_ca
                                # construct the first term using the first copy
                                # replace P2_r in the copy list with P2_{index from left P2 in commutator}
                                sub_terms_list_copy_1[sub_terms_list_copy_1.index(P2_r)] = f"P2{left_P2_commutator_index}"
                                # then delete the P2 that was in the LHS of the commutator the copy list since this is accounted for in the commutator
                                left_P2_commutator_pos = sub_terms_list_copy_1.index(P2_var_on_right_commute_list)
                                sub_terms_list_copy_1.pop(left_P2_commutator_pos)
                                # then also delete the P2 that was in the RHS of the commutator from the list
                                right_P2_commutator_pos = sub_terms_list_copy_1.index(P_var)
                                sub_terms_list_copy_1.pop(right_P2_commutator_pos)
                                # then in this term, we will also have a metric to contract: g_{P2_r index}_{index from right P2 in commutator}
                                # so find another term with one of these indices (choose P2_r cos why not) and change its letter index accordingly
                                if P2_var_on_right_commute_list[-1] == P2_r[-1]:
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['d*(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))
                                else:
                                    for var in sub_terms_list_copy_1:
                                        # check if index is same as P2_r -> this var will be contracted with the metric
                                        if var[-1] == P2_r[-1]:
                                            # find the index of this var in the list copy
                                            sub_list_var_index = sub_terms_list_copy_1.index(var)
                                            # if the P2_r index is being contracted with var, then the remaining index in the index of the
                                            # P2 in the RHS of the commutator
                                            sub_terms_list_copy_1[sub_list_var_index] = f"{var[:-2]}{right_P2_commutator_index}"
                                            break
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))

                                # now repeat the process but for the second term; same structure, except for a sign and swap of indices
                                # construct the second term using the second copy
                                # delete the P2_r in the copy list since this is accounted for in the commutator
                                left_P2_commutator_pos = sub_terms_list_copy_2.index(P2_var_on_right_commute_list)
                                sub_terms_list_copy_2.pop(left_P2_commutator_pos)
                                # then also delete the P2 that was in the RHS of the commutator from the list
                                right_P2_commutator_pos = sub_terms_list_copy_2.index(P_var)
                                sub_terms_list_copy_2.pop(right_P2_commutator_pos)
                                # replace P2_r in the copy list with P2_{index from right P2 in commutator}
                                sub_terms_list_copy_2[sub_terms_list_copy_2.index(P2_r)] = f"P2{right_P2_commutator_index}"
                                # then in this term, we will also have a metric to contract: g_{P2_r index}_{index from left P2 in commutator}
                                # so find another term with one of these indices (choose P2_r cos why not) and change its letter index accordingly
                                if P_var[-1] == P2_r[-1]:
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-d*(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                                else:
                                    for var in sub_terms_list_copy_2:
                                        # check if index is same as P2_r -> this var will be contracted with the metric
                                        if var[-1] == P2_r[-1]:
                                            # find the index of this var in the list copy
                                            sub_list_var_index = sub_terms_list_copy_2.index(var)
                                            # if the P2_r index is being contracted with var, then the remaining index in the index of the
                                            # P2 in the LHS of the commutator
                                            sub_terms_list_copy_2[sub_list_var_index] = f"{var[:-2]}{left_P2_commutator_index}"
                                            break
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                            
                            # now we have to add in the riemann tensor operator that acts on the fields themselves.
                            # again, make two separate terms
                            sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*a2_b*U2_a
                            sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*a2_a*U2_b
                            # the indices are now only taken from P_var and P_var_on_left
                            # replace the LHS commutator P1 term with aU term
                            LHS_commutator_P2_index = sub_terms_list_copy_1.index(P2_var_on_right_commute_list)
                            sub_terms_list_copy_1[LHS_commutator_P2_index] = f"a2{left_P2_commutator_index}*U2{right_P2_commutator_index}"
                            # then delete the RHS commutator P2 term from the list
                            RHS_commutator_P2_index = sub_terms_list_copy_1.index(P_var)
                            sub_terms_list_copy_1.pop(RHS_commutator_P2_index)
                            # then add in the 1/l^2 term, accounting for signs
                            if sub_terms_list_copy_1[0][0] != '-':
                                sub_terms_list_copy_1 = ['-(1/l^2)'] + sub_terms_list_copy_1
                            else:
                                sub_terms_list_copy_1 = ['(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_1))

                            # then the second term
                            # replace the LHS commutator P2 term with aU term
                            LHS_commutator_P2_index = sub_terms_list_copy_2.index(P2_var_on_right_commute_list)
                            sub_terms_list_copy_2[LHS_commutator_P2_index] = f"a2{right_P2_commutator_index}*U2{left_P2_commutator_index}"
                            # then delete the RHS commutator P2 term from the list
                            RHS_commutator_P2_index = sub_terms_list_copy_2.index(P_var)
                            sub_terms_list_copy_2.pop(RHS_commutator_P2_index)
                            # then add in the 1/l^2 term, accounting for signs
                            if sub_terms_list_copy_2[0][0] != '-':
                                sub_terms_list_copy_2 = ['(1/l^2)'] + sub_terms_list_copy_2
                            else:
                                sub_terms_list_copy_2 = ['-(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # lastly swap the orders of P1s in both lists for the next iteration
                            sub_terms_list[P2_var_on_right_sub_list], sub_terms_list[P2_index_sub_list] = sub_terms_list[P2_index_sub_list], sub_terms_list[P2_var_on_right_sub_list]
                            P2_to_commute[P2_index_commute_list+1], P2_to_commute[P2_index_commute_list] = P2_to_commute[P2_index_commute_list], P2_to_commute[P2_index_commute_list+1]    
                        # only want to do this for the leftmost U1P3
                        break
        if not U2P2_present:
            # just reconstruct the original term
            expanded_terms.append('*'.join(sub_terms_list))
    
    return pull_minus_signs_to_front(' + '.join(expanded_terms))

def pull_Div3_right(equation):
    # pulls the leftmost U3P3 to the right to discard it
    # break the equation down into individual terms by top level '+' symbols
    equation = replace_minuses(equation)
    terms = top_level_split(equation)
    expanded_terms = []
    for term in terms:
        # break into individual variables by * signs
        # this will be the main list which will be rearranged as we go
        sub_terms_list = term.split('*')
        # check over sub_terms_list for U3 variables
        U3P3_present = False
        for U_var in reversed(sub_terms_list):
            if U3P3_present:
                break
            elif U_var.startswith('U3'):
                # store the letter index and position index of the U3
                U3_lett_index = U_var[-1] 
                U3_pos_index = sub_terms_list.index(U_var)
                for P_var in sub_terms_list:
                    if P_var.startswith('P3') and P_var[-1] == U3_lett_index:
                        U3P3_present = True
                        # store the P3 index (letter index should be same as U3 if relevant)
                        P3_pos_index = sub_terms_list.index(P_var)
                        # we have found a contracted U3P3, so then we can do the rearrangement
                        # now we want to move all the U3, Z1, Z2 left of the leftmost P3 before we can commute the P3s
                        # get all the U3, Z1, Z2 variables in a list
                        U3_Z1_Z2_list = [var for var in sub_terms_list if var.startswith('U1') or var.startswith('U2') or var.startswith('U3') or var.startswith('Z2') or var.startswith('Z1') or var.startswith('Z3')]
                        # get the list without these variables
                        remaining_list = [var for var in sub_terms_list if var not in U3_Z1_Z2_list]
                        # get the position of the leftmost P3 in the remaining list
                        for var in remaining_list:
                            if var.startswith('P3'):
                                leftmost_P3 = var
                                leftmost_P3_pos = remaining_list.index(var)
                                break
                        # now sequentially insert the U2_Z1_Z3_list on the left of the leftmost P3
                        for U3_Z1_Z2 in U3_Z1_Z2_list:
                            remaining_list.insert(leftmost_P3_pos, U3_Z1_Z2)
                        # assign rearranged list to our main sub_terms_list
                        sub_terms_list = remaining_list
                        # now we can start commuting P3s. We make a list of all the P3s
                        P3_to_commute = [var for var in sub_terms_list if var.startswith('P3')]
                        # want our P3 in the U3P3 contraction to be moved right - if not in final position of P3_to_commute, carry on
                        # rearranging 
                        while P3_to_commute.index(P_var) != (len(P3_to_commute)-1):
                            # find P3 index in the to_commute list
                            P3_index_commute_list = P3_to_commute.index(P_var)
                            # find P3 index in the sub_terms_list
                            P3_index_sub_list = sub_terms_list.index(P_var)
                            # store the P3 variable on the next right to be commuted through
                            P3_var_on_right_commute_list = P3_to_commute[P3_index_commute_list+1]
                            # find this variable's position in the sub_terms_list
                            P3_var_on_right_sub_list = sub_terms_list.index(P3_var_on_right_commute_list)
                            # extract relevant information about the P3s being commuted (indices)
                            # used for constructing the indices in the riemann tensor
                            right_P3_commutator_index = P3_var_on_right_commute_list[-2:]
                            left_P3_commutator_index = P_var[-2:]
                            # create terms arising from the riemann tensor contracting with P3s on the right
                            P3s_on_right = [var for var in P3_to_commute if P3_to_commute.index(var) > P3_index_commute_list+1]
                            for P3_r in P3s_on_right:
                                # create copy to manipulate terms without altering master list
                                sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*P3_a*g_cb
                                sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*P3_b*g_ca
                                # construct the first term using the first copy
                                # replace P3_r in the copy list with P3_{index from left P3 in commutator}
                                sub_terms_list_copy_1[sub_terms_list_copy_1.index(P3_r)] = f"P3{left_P3_commutator_index}"
                                # then delete the P3 that was in the LHS of the commutator the copy list since this is accounted for in the commutator
                                left_P3_commutator_pos = sub_terms_list_copy_1.index(P3_var_on_right_commute_list)
                                sub_terms_list_copy_1.pop(left_P3_commutator_pos)
                                # then also delete the P3 that was in the RHS of the commutator from the list
                                right_P3_commutator_pos = sub_terms_list_copy_1.index(P_var)
                                sub_terms_list_copy_1.pop(right_P3_commutator_pos)
                                # then in this term, we will also have a metric to contract: g_{P3_r index}_{index from right P3 in commutator}
                                # so find another term with one of these indices (choose P3_r cos why not) and change its letter index accordingly
                                if P3_var_on_right_commute_list[-1] == P3_r[-1]:
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['d*(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))
                                else:
                                    for var in sub_terms_list_copy_1:
                                        # check if index is same as P3_r -> this var will be contracted with the metric
                                        if var[-1] == P3_r[-1]:
                                            # find the index of this var in the list copy
                                            sub_list_var_index = sub_terms_list_copy_1.index(var)
                                            # if the P3_r index is being contracted with var, then the remaining index in the index of the
                                            # P3 in the RHS of the commutator
                                            sub_terms_list_copy_1[sub_list_var_index] = f"{var[:-2]}{right_P3_commutator_index}"
                                            break
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))

                                # now repeat the process but for the second term; same structure, except for a sign and swap of indices
                                # construct the second term using the second copy
                                # delete the P3_r in the copy list since this is accounted for in the commutator
                                left_P3_commutator_pos = sub_terms_list_copy_2.index(P3_var_on_right_commute_list)
                                sub_terms_list_copy_2.pop(left_P3_commutator_pos)
                                # then also delete the P3 that was in the RHS of the commutator from the list
                                right_P3_commutator_pos = sub_terms_list_copy_2.index(P_var)
                                sub_terms_list_copy_2.pop(right_P3_commutator_pos)
                                # replace P3_r in the copy list with P3_{index from right P3 in commutator}
                                sub_terms_list_copy_2[sub_terms_list_copy_2.index(P3_r)] = f"P3{right_P3_commutator_index}"
                                # then in this term, we will also have a metric to contract: g_{P3_r index}_{index from left P3 in commutator}
                                # so find another term with one of these indices (choose P3_r cos why not) and change its letter index accordingly
                                if P_var[-1] == P3_r[-1]:
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-d*(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                                else:
                                    for var in sub_terms_list_copy_2:
                                        # check if index is same as P3_r -> this var will be contracted with the metric
                                        if var[-1] == P3_r[-1]:
                                            # find the index of this var in the list copy
                                            sub_list_var_index = sub_terms_list_copy_2.index(var)
                                            # if the P3_r index is being contracted with var, then the remaining index in the index of the
                                            # P3 in the LHS of the commutator
                                            sub_terms_list_copy_2[sub_list_var_index] = f"{var[:-2]}{left_P3_commutator_index}"
                                            break
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                            
                            # now we have to add in the riemann tensor operator that acts on the fields themselves.
                            # again, make two separate terms
                            sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*a3_b*U3_a
                            sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*a3_a*U3_b
                            # the indices are now only taken from P_var and P_var_on_left
                            # replace the LHS commutator P1 term with aU term
                            LHS_commutator_P3_index = sub_terms_list_copy_1.index(P3_var_on_right_commute_list)
                            sub_terms_list_copy_1[LHS_commutator_P3_index] = f"a3{left_P3_commutator_index}*U3{right_P3_commutator_index}"
                            # then delete the RHS commutator P3 term from the list
                            RHS_commutator_P3_index = sub_terms_list_copy_1.index(P_var)
                            sub_terms_list_copy_1.pop(RHS_commutator_P3_index)
                            # then add in the 1/l^2 term, accounting for signs
                            if sub_terms_list_copy_1[0][0] != '-':
                                sub_terms_list_copy_1 = ['-(1/l^2)'] + sub_terms_list_copy_1
                            else:
                                sub_terms_list_copy_1 = ['(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_1))

                            # then the second term
                            # replace the LHS commutator P3 term with aU term
                            LHS_commutator_P3_index = sub_terms_list_copy_2.index(P3_var_on_right_commute_list)
                            sub_terms_list_copy_2[LHS_commutator_P3_index] = f"a3{right_P3_commutator_index}*U3{left_P3_commutator_index}"
                            # then delete the RHS commutator P3 term from the list
                            RHS_commutator_P3_index = sub_terms_list_copy_2.index(P_var)
                            sub_terms_list_copy_2.pop(RHS_commutator_P3_index)
                            # then add in the 1/l^2 term, accounting for signs
                            if sub_terms_list_copy_2[0][0] != '-':
                                sub_terms_list_copy_2 = ['(1/l^2)'] + sub_terms_list_copy_2
                            else:
                                sub_terms_list_copy_2 = ['-(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # lastly swap the orders of P1s in both lists for the next iteration
                            sub_terms_list[P3_var_on_right_sub_list], sub_terms_list[P3_index_sub_list] = sub_terms_list[P3_index_sub_list], sub_terms_list[P3_var_on_right_sub_list]
                            P3_to_commute[P3_index_commute_list+1], P3_to_commute[P3_index_commute_list] = P3_to_commute[P3_index_commute_list], P3_to_commute[P3_index_commute_list+1]    
                        # only want to do this for the leftmost U3P3
                        break
        if not U3P3_present:
            # just reconstruct the original term
            expanded_terms.append('*'.join(sub_terms_list))
    
    return pull_minus_signs_to_front(' + '.join(expanded_terms))

def remove_traces(equation):
    """
    Removes trace terms from an equation string. A trace term is one where there are contracted U variables.
    
    Parameters:
        equation (str): The input equation as a single string.

    Returns:
        str: The equation with trace terms removed.
    """
    equation = reorder_un_terms(equation)
    # Split the equation into individual terms
    terms = top_level_split(equation)
    # Initialize a list to store valid terms
    valid_terms = []
    
    # Regular expressions for detecting contracted U variables
    contraction_patterns = [
        r'U(\d)_([a-z])\*U\1\^\2',  # Pattern for contracted subscripts/superscripts
        r'U(\d)\^([a-z])\*U\1_\2'   # Swapped order
    ]
    
    # Process each term
    for term in terms:
        is_trace = False
        # Check all contraction patterns
        for pattern in contraction_patterns:
            if re.search(pattern, term):
                is_trace = True
                break  # Skip this term if it's a trace term
        if not is_trace:
            valid_terms.append(term)
    
    # Recombine valid terms back into the original form
    reconstructed_equation = ' + '.join(valid_terms)
    return reconstructed_equation

def pull_all_Divs_right(equation):
    return move_all_a_left(pull_Div1_right(move_all_a_left(pull_Div2_right(move_all_a_left(pull_Div3_right(move_all_a_left(equation)))))))

def pull_P3P1s_left(equation):
    # pulls the leftmost P3P1 to the left to integrate by parts
    # break the equation down into individual terms by top level '+' symbols
    equation = replace_minuses(equation)
    terms = top_level_split(equation)
    expanded_terms = []
    for term in terms:
        # break into individual variables by * signs
        # this will be the main list which will be rearranged as we go
        sub_terms_list = term.split('*')
        # check over sub_terms_list for P3 variables
        P3P1_present = False
        for P3_var in sub_terms_list:
            if P3P1_present:
                break
            elif P3_var.startswith('P3'):
                P3_var_stored = P3_var
                # store the letter index and position index of the P3
                P3_lett_index = P3_var[-1] 
                P3_pos_index = sub_terms_list.index(P3_var)
                for P1_var in sub_terms_list:
                    if P1_var.startswith('P1') and P1_var[-1] == P3_lett_index:
                        P3P1_present = True
                        # store the P1 index (letter index should be same as P3 if relevant)
                        P1_pos_index = sub_terms_list.index(P1_var)
                        # we have found a contracted P3P1, so then we can do the rearrangement
                        # now we want to move all the U1, U3, Z1, Z2, Z3 left of the leftmost P before we can commute P1s/P3s
                        # get all the U1, U3, Z1, Z2, Z3 variables in a list
                        U1_U3_Z1_Z2_Z3_list = [var for var in sub_terms_list if var.startswith('U1') or var.startswith('U2') or var.startswith('U3') or var.startswith('Z2') or var.startswith('Z1') or var.startswith('Z3')]
                        # get the list without these variables
                        remaining_list = [var for var in sub_terms_list if var not in U1_U3_Z1_Z2_Z3_list]
                        # get the position of the leftmost P in the remaining list
                        for var in remaining_list:
                            if var.startswith('P'):
                                leftmost_P = var
                                leftmost_P_pos = remaining_list.index(var)
                                break
                        # now sequentially insert the U1_U3_Z1_Z2_Z3_list on the left of the leftmost P
                        for U1_U3_Z1_Z2_Z3 in U1_U3_Z1_Z2_Z3_list:
                            remaining_list.insert(leftmost_P_pos, U1_U3_Z1_Z2_Z3)
                        # assign rearranged list to our main sub_terms_list
                        sub_terms_list = remaining_list
                        # now we can start commuting P3s. We make a list of all the P3s
                        P3_to_commute = [var for var in sub_terms_list if var.startswith('P3')]
                        # want our P3 in the U1P3 contraction to be moved left - if not in position zero of P3_to_commute, carry on
                        # rearranging
                        while P3_to_commute.index(P3_var) != 0:
                            # find P3 index in the to_commute list
                            P3_index_commute_list = P3_to_commute.index(P3_var)
                            # find P3 index in the sub_terms_list
                            P3_index_sub_list = sub_terms_list.index(P3_var)
                            # store the P3 variable on the next left to be commuted through
                            P3_var_on_left_commute_list = P3_to_commute[P3_index_commute_list-1]
                            # find this variable's position in the sub_terms_list
                            P3_var_on_left_sub_list = sub_terms_list.index(P3_var_on_left_commute_list)
                            # extract relevant information about the P3s being commuted (indices)
                            # used for constructing the indices in the riemann tensor
                            left_P3_commutator_index = P3_var_on_left_commute_list[-2:]
                            right_P3_commutator_index = P3_var[-2:]
                            # create terms arising from the riemann tensor contracting with P3s on the right
                            P3s_on_right = [var for var in P3_to_commute if P3_to_commute.index(var) > P3_index_commute_list]
                            for P3_r in P3s_on_right:
                                # create copy to manipulate terms without altering master list
                                sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*P1_a*g_cb
                                sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*P1_b*g_ca
                                # construct the first term using the first copy
                                # replace P3_r in the copy list with P3_{index from left P3 in commutator}
                                sub_terms_list_copy_1[sub_terms_list_copy_1.index(P3_r)] = f"P3{left_P3_commutator_index}"
                                # then delete the P3 that was in the LHS of the commutator the copy list since this is accounted for in the commutator
                                left_P3_commutator_pos = sub_terms_list_copy_1.index(P3_var_on_left_commute_list)
                                sub_terms_list_copy_1.pop(left_P3_commutator_pos)
                                # then also delete the P3 that was in the RHS of the commutator from the list
                                right_P3_commutator_pos = sub_terms_list_copy_1.index(P3_var)
                                sub_terms_list_copy_1.pop(right_P3_commutator_pos)
                                # then in this term, we will also have a metric to contract: g_{P3_r index}_{index from right P3 in commutator}
                                # so find another term with one of these indices (choose P3_r cos why not) and change its letter index accordingly
                                if P3_var[-1] == P3_r[-1]:
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['d*(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))
                                else:
                                    for var in sub_terms_list_copy_1:
                                        # check if index is same as P3_r -> this var will be contracted with the metric
                                        if var[-1] == P3_r[-1]:
                                            # find the index of this var in the list copy
                                            sub_list_var_index = sub_terms_list_copy_1.index(var)
                                            # if the P3_r index is being contracted with var, then the remaining index in the index of the
                                            # P3 in the RHS of the commutator
                                            sub_terms_list_copy_1[sub_list_var_index] = f"{var[:-2]}{right_P3_commutator_index}"
                                            break
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))

                                # now repeat the process but for the second term; same structure, except for a sign and swap of indices
                                # construct the second term using the second copy
                                # delete the P3_r in the copy list since this is accounted for in the commutator
                                left_P3_commutator_pos = sub_terms_list_copy_2.index(P3_var_on_left_commute_list)
                                sub_terms_list_copy_2.pop(left_P3_commutator_pos)
                                # then also delete the P3 that was in the RHS of the commutator from the list
                                right_P3_commutator_pos = sub_terms_list_copy_2.index(P3_var)
                                sub_terms_list_copy_2.pop(right_P3_commutator_pos)
                                # replace P3_r in the copy list with P3_{index from right P3 in commutator}
                                sub_terms_list_copy_2[sub_terms_list_copy_2.index(P3_r)] = f"P3{right_P3_commutator_index}"
                                # then in this term, we will also have a metric to contract: g_{P3_r index}_{index from left P3 in commutator}
                                # so find another term with one of these indices (choose P3_r cos why not) and change its letter index accordingly
                                if P3_var_on_left_commute_list[-1] == P3_r[-1]:
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-d*(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                                else:
                                    for var in sub_terms_list_copy_2:
                                        # check if index is same as P3_r -> this var will be contracted with the metric
                                        if var[-1] == P3_r[-1]:
                                            # find the index of this var in the list copy
                                            sub_list_var_index = sub_terms_list_copy_2.index(var)
                                            # if the P3_r index is being contracted with var, then the remaining index in the index of the
                                            # P3 in the LHS of the commutator
                                            sub_terms_list_copy_2[sub_list_var_index] = f"{var[:-2]}{left_P3_commutator_index}"
                                            break
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                            
                            # now we have to add in the riemann tensor operator that acts on the fields themselves.
                            # again, make two separate terms
                            sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*a3_b*U3_a
                            sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*a3_a*U3_b
                            # the indices are now only taken from P_var and P_var_on_left
                            # replace the LHS commutator P3 term with aU term
                            LHS_commutator_P3_index = sub_terms_list_copy_1.index(P3_var_on_left_commute_list)
                            sub_terms_list_copy_1[LHS_commutator_P3_index] = f"a3{left_P3_commutator_index}*U3{right_P3_commutator_index}"
                            # then delete the RHS commutator P3 term from the list
                            RHS_commutator_P3_index = sub_terms_list_copy_1.index(P3_var)
                            sub_terms_list_copy_1.pop(RHS_commutator_P3_index)
                            # then add in the 1/l^2 term, accounting for signs
                            if sub_terms_list_copy_1[0][0] != '-':
                                sub_terms_list_copy_1 = ['-(1/l^2)'] + sub_terms_list_copy_1
                            else:
                                sub_terms_list_copy_1 = ['(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_1))

                            # then the second term
                            # replace the LHS commutator P3 term with aU term
                            LHS_commutator_P3_index = sub_terms_list_copy_2.index(P3_var_on_left_commute_list)
                            sub_terms_list_copy_2[LHS_commutator_P3_index] = f"a3{right_P3_commutator_index}*U3{left_P3_commutator_index}"
                            # then delete the RHS commutator P3 term from the list
                            RHS_commutator_P3_index = sub_terms_list_copy_2.index(P3_var)
                            sub_terms_list_copy_2.pop(RHS_commutator_P3_index)
                            # then add in the 1/l^2 term, accounting for signs
                            if sub_terms_list_copy_2[0][0] != '-':
                                sub_terms_list_copy_2 = ['(1/l^2)'] + sub_terms_list_copy_2
                            else:
                                sub_terms_list_copy_2 = ['-(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # lastly swap the orders of P3s in both lists for the next iteration
                            sub_terms_list[P3_var_on_left_sub_list], sub_terms_list[P3_index_sub_list] = sub_terms_list[P3_index_sub_list], sub_terms_list[P3_var_on_left_sub_list]
                            P3_to_commute[P3_index_commute_list-1], P3_to_commute[P3_index_commute_list] = P3_to_commute[P3_index_commute_list], P3_to_commute[P3_index_commute_list-1]

                        # now we can start commuting P1s. We make a list of all the P1s
                        P1_to_commute = [var for var in sub_terms_list if var.startswith('P1')]
                        # want our P1 in the P3P1 contraction to be moved left - if not in position zero of P1_to_commute, carry on
                        # rearranging 
                        while P1_to_commute.index(P1_var) != 0:
                            # find P1 index in the to_commute list
                            P1_index_commute_list = P1_to_commute.index(P1_var)
                            # find P1 index in the sub_terms_list
                            P1_index_sub_list = sub_terms_list.index(P1_var)
                            # store the P1 variable on the next left to be commuted through
                            P1_var_on_left_commute_list = P1_to_commute[P1_index_commute_list-1]
                            # find this variable's position in the sub_terms_list
                            P1_var_on_left_sub_list = sub_terms_list.index(P1_var_on_left_commute_list)
                            # extract relevant information about the P1s being commuted (indices)
                            # used for constructing the indices in the riemann tensor
                            left_P1_commutator_index = P1_var_on_left_commute_list[-2:]
                            right_P1_commutator_index = P1_var[-2:]
                            # create terms arising from the riemann tensor contracting with P1s on the right
                            P1s_on_right = [var for var in P1_to_commute if P1_to_commute.index(var) > P1_index_commute_list]
                            for P1_r in P1s_on_right:
                                # create copy to manipulate terms without altering master list
                                sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*P1_a*g_cb
                                sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*P1_b*g_ca
                                # construct the first term using the first copy
                                # replace P3_r in the copy list with P1_{index from left P1 in commutator}
                                sub_terms_list_copy_1[sub_terms_list_copy_1.index(P1_r)] = f"P1{left_P1_commutator_index}"
                                # then delete the P1 that was in the LHS of the commutator the copy list since this is accounted for in the commutator
                                left_P1_commutator_pos = sub_terms_list_copy_1.index(P1_var_on_left_commute_list)
                                sub_terms_list_copy_1.pop(left_P1_commutator_pos)
                                # then also delete the P1 that was in the RHS of the commutator from the list
                                right_P1_commutator_pos = sub_terms_list_copy_1.index(P1_var)
                                sub_terms_list_copy_1.pop(right_P1_commutator_pos)
                                # then in this term, we will also have a metric to contract: g_{P1_r index}_{index from right P1 in commutator}
                                # so find another term with one of these indices (choose P1_r cos why not) and change its letter index accordingly
                                if P1_var[-1] == P1_r[-1]:
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['d*(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))
                                else:
                                    for var in sub_terms_list_copy_1:
                                        # check if index is same as P1_r -> this var will be contracted with the metric
                                        if var[-1] == P1_r[-1]:
                                            # find the index of this var in the list copy
                                            sub_list_var_index = sub_terms_list_copy_1.index(var)
                                            # if the P1_r index is being contracted with var, then the remaining index in the index of the
                                            # P1 in the RHS of the commutator
                                            sub_terms_list_copy_1[sub_list_var_index] = f"{var[:-2]}{right_P1_commutator_index}"
                                            break
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))

                                # now repeat the process but for the second term; same structure, except for a sign and swap of indices
                                # construct the second term using the second copy
                                # delete the P1_r in the copy list since this is accounted for in the commutator
                                left_P1_commutator_pos = sub_terms_list_copy_2.index(P1_var_on_left_commute_list)
                                sub_terms_list_copy_2.pop(left_P1_commutator_pos)
                                # then also delete the P1 that was in the RHS of the commutator from the list
                                right_P1_commutator_pos = sub_terms_list_copy_2.index(P1_var)
                                sub_terms_list_copy_2.pop(right_P1_commutator_pos)
                                # replace P1_r in the copy list with P1_{index from right P1 in commutator}
                                sub_terms_list_copy_2[sub_terms_list_copy_2.index(P1_r)] = f"P1{right_P1_commutator_index}"
                                # then in this term, we will also have a metric to contract: g_{P1_r index}_{index from left P1 in commutator}
                                # so find another term with one of these indices (choose P1_r cos why not) and change its letter index accordingly
                                if P1_var_on_left_commute_list[-1] == P1_r[-1]:
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-d*(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                                else:
                                    for var in sub_terms_list_copy_2:
                                        # check if index is same as P1_r -> this var will be contracted with the metric
                                        if var[-1] == P1_r[-1]:
                                            # find the index of this var in the list copy
                                            sub_list_var_index = sub_terms_list_copy_2.index(var)
                                            # if the P1_r index is being contracted with var, then the remaining index in the index of the
                                            # P1 in the LHS of the commutator
                                            sub_terms_list_copy_2[sub_list_var_index] = f"{var[:-2]}{left_P1_commutator_index}"
                                            break
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                            
                            # now we have to add in the riemann tensor operator that acts on the fields themselves.
                            # again, make two separate terms
                            sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*a1_b*U1_a
                            sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*a1_a*U1_b
                            # the indices are now only taken from P_var and P_var_on_left
                            # replace the LHS commutator P1 term with aU term
                            LHS_commutator_P1_index = sub_terms_list_copy_1.index(P1_var_on_left_commute_list)
                            sub_terms_list_copy_1[LHS_commutator_P1_index] = f"a1{left_P1_commutator_index}*U1{right_P1_commutator_index}"
                            # then delete the RHS commutator P1 term from the list
                            RHS_commutator_P1_index = sub_terms_list_copy_1.index(P1_var)
                            sub_terms_list_copy_1.pop(RHS_commutator_P1_index)
                            # then add in the 1/l^2 term, accounting for signs
                            if sub_terms_list_copy_1[0][0] != '-':
                                sub_terms_list_copy_1 = ['-(1/l^2)'] + sub_terms_list_copy_1
                            else:
                                sub_terms_list_copy_1 = ['(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_1))

                            # then the second term
                            # replace the LHS commutator P1 term with aU term
                            LHS_commutator_P1_index = sub_terms_list_copy_2.index(P1_var_on_left_commute_list)
                            sub_terms_list_copy_2[LHS_commutator_P1_index] = f"a1{right_P1_commutator_index}*U1{left_P1_commutator_index}"
                            # then delete the RHS commutator P1 term from the list
                            RHS_commutator_P1_index = sub_terms_list_copy_2.index(P1_var)
                            sub_terms_list_copy_2.pop(RHS_commutator_P1_index)
                            # then add in the 1/l^2 term, accounting for signs
                            if sub_terms_list_copy_2[0][0] != '-':
                                sub_terms_list_copy_2 = ['(1/l^2)'] + sub_terms_list_copy_2
                            else:
                                sub_terms_list_copy_2 = ['-(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # lastly swap the orders of P3s in both lists for the next iteration
                            sub_terms_list[P1_var_on_left_sub_list], sub_terms_list[P1_index_sub_list] = sub_terms_list[P1_index_sub_list], sub_terms_list[P1_var_on_left_sub_list]
                            P1_to_commute[P1_index_commute_list-1], P1_to_commute[P1_index_commute_list] = P1_to_commute[P1_index_commute_list], P1_to_commute[P1_index_commute_list-1]
                        # only want to do this for the leftmost P3P1
                        break
        if P3P1_present:
            # now that sub_terms_list has the variables in the desired order, we can replace P3P1 with (1/2)*(P2_a*P2^a - P1_a*P1^a - P3_a*P3^a)
            # create copies again to construct desired terms
            sub_terms_list_B1 = sub_terms_list.copy()
            sub_terms_list_B2 = sub_terms_list.copy()
            sub_terms_list_B3 = sub_terms_list.copy()
            # can just use P3_var and P1_var - in first list, delete the P3_var and replace the P1_var with the box operator (same index)
            # but need the replaced terms to be on the left of all derivatives due to the partial integration
            P1_var_index = sub_terms_list.index(P1_var)
            leftmost_P_index = P1_var_index
            for var in sub_terms_list_B1:
                # find left most P index
                if var.startswith('P'):
                    leftmost_P_index = sub_terms_list_B1.index(var)
                    break
            # then replace with relevant box operators
            lett_index_to_use = P3_var_stored[-1]
            sub_terms_list_B1.insert(leftmost_P_index, f'P1_{lett_index_to_use}*P1^{lett_index_to_use}')
            sub_terms_list_B2.insert(leftmost_P_index, f'P2_{lett_index_to_use}*P2^{lett_index_to_use}')
            sub_terms_list_B3.insert(leftmost_P_index, f'P3_{lett_index_to_use}*P3^{lett_index_to_use}')
            # B2 term is multiplied by 1/2, while B1 and B3 are multtiplied by -1/2 

            P3_var_index_B1 = sub_terms_list_B1.index(P3_var_stored)
            P3_var_index_B2 = sub_terms_list_B2.index(P3_var_stored)
            P3_var_index_B3 = sub_terms_list_B3.index(P3_var_stored)
            # delete P3_var in all lists
            sub_terms_list_B1.pop(P3_var_index_B1)
            sub_terms_list_B2.pop(P3_var_index_B2)
            sub_terms_list_B3.pop(P3_var_index_B3)

            # and delete the P1_var too
            P1_var_index_B1 = sub_terms_list_B1.index(P1_var)
            P1_var_index_B2 = sub_terms_list_B2.index(P1_var)
            P1_var_index_B3 = sub_terms_list_B3.index(P1_var)
            sub_terms_list_B1.pop(P1_var_index_B1)
            sub_terms_list_B2.pop(P1_var_index_B2)
            sub_terms_list_B3.pop(P1_var_index_B3)

            if sub_terms_list_B1[0][0] == '-':
                sub_terms_list_B1 = ['(1/2)'] + [sub_terms_list_B1[0][1:]] + sub_terms_list_B1[1:]
            else:
                sub_terms_list_B1 = ['-(1/2)'] + sub_terms_list_B1
            if sub_terms_list_B2[0][0] == '-':
                sub_terms_list_B2 = ['-(1/2)'] + [sub_terms_list_B2[0][1:]] + sub_terms_list_B2[1:]
            else:
                sub_terms_list_B2 = ['(1/2)'] + sub_terms_list_B2
            if sub_terms_list_B3[0][0] == '-':
                sub_terms_list_B3 = ['(1/2)'] + [sub_terms_list_B3[0][1:]] + sub_terms_list_B3[1:]
            else:
                sub_terms_list_B3 = ['-(1/2)'] + sub_terms_list_B3
            
            # then we can construct the terms by joining with *
            expanded_terms.append('*'.join(sub_terms_list_B1))
            expanded_terms.append('*'.join(sub_terms_list_B2))
            expanded_terms.append('*'.join(sub_terms_list_B3))
        else:
            # just reconstruct the original term
            expanded_terms.append('*'.join(sub_terms_list))

    return pull_minus_signs_to_front(' + '.join(expanded_terms))

def pull_P3P2s_left(equation):
    # pulls the leftmost P3P2 to the left to integrate by parts
    # break the equation down into individual terms by top level '+' symbols
    equation = replace_minuses(equation)
    terms = top_level_split(equation)
    expanded_terms = []
    for term in terms:
        # break into individual variables by * signs
        # this will be the main list which will be rearranged as we go
        sub_terms_list = term.split('*')
        # check over sub_terms_list for P3 variables
        P3P2_present = False
        for P3_var in sub_terms_list:
            if P3P2_present == True:
                break
            elif P3_var.startswith('P3'):
                P3_var_stored = P3_var
                # store the letter index and position index of the P3
                P3_lett_index = P3_var[-1]
                P3_pos_index = sub_terms_list.index(P3_var)
                for P2_var in sub_terms_list:
                    if P2_var.startswith('P2') and P2_var[-1] == P3_lett_index:
                        P3P2_present = True
                        # store the P2 index (letter index should be same as P3 if relevant)
                        P2_pos_index = sub_terms_list.index(P2_var)

                        # -------------------------------------------------------
                        # Now we want to move all U2, U3, Z1, Z2, Z3 left of the leftmost P
                        # before we can commute P2s/P3s
                        # -------------------------------------------------------
                        U2_U3_Z1_Z2_Z3_list = [var for var in sub_terms_list if var.startswith('U1') or var.startswith('U2') or var.startswith('U3') or var.startswith('Z2') or var.startswith('Z1') or var.startswith('Z3')]
                        # Create a list with those items removed
                        remaining_list = [var for var in sub_terms_list if var not in U2_U3_Z1_Z2_Z3_list]

                        # Find the position of the leftmost P in the remaining_list
                        for var in remaining_list:
                            if var.startswith('P'):
                                leftmost_P = var
                                leftmost_P_pos = remaining_list.index(var)
                                break

                        # Insert all U2/U3/Z1/Z2/Z3 vars to the left of that leftmost P
                        for U2_U3_Z1_Z2_Z3 in U2_U3_Z1_Z2_Z3_list:
                            remaining_list.insert(leftmost_P_pos, U2_U3_Z1_Z2_Z3)

                        # Reassign sub_terms_list
                        sub_terms_list = remaining_list

                        # -------------------------------------------------------
                        # Commute P3s
                        # -------------------------------------------------------
                        P3_to_commute = [var for var in sub_terms_list if var.startswith('P3')]
                        while P3_to_commute.index(P3_var) != 0:
                            P3_index_commute_list = P3_to_commute.index(P3_var)
                            P3_index_sub_list = sub_terms_list.index(P3_var)
                            P3_var_on_left_commute_list = P3_to_commute[P3_index_commute_list - 1]
                            P3_var_on_left_sub_list = sub_terms_list.index(P3_var_on_left_commute_list)
                            left_P3_commutator_index = P3_var_on_left_commute_list[-2:]
                            right_P3_commutator_index = P3_var[-2:]

                            # Create terms from Riemann tensor contracting with P3s on the right
                            P3s_on_right = [
                                var for var in P3_to_commute
                                if P3_to_commute.index(var) > P3_index_commute_list
                            ]
                            for P3_r in P3s_on_right:
                                # copy 1
                                sub_terms_list_copy_1 = sub_terms_list.copy()
                                # copy 2
                                sub_terms_list_copy_2 = sub_terms_list.copy()

                                # First term
                                sub_terms_list_copy_1[sub_terms_list_copy_1.index(P3_r)] = f"P3{left_P3_commutator_index}"
                                left_P3_commutator_pos = sub_terms_list_copy_1.index(P3_var_on_left_commute_list)
                                sub_terms_list_copy_1.pop(left_P3_commutator_pos)
                                right_P3_commutator_pos = sub_terms_list_copy_1.index(P3_var)
                                sub_terms_list_copy_1.pop(right_P3_commutator_pos)
                                if P3_var[-1] == P3_r[-1]:
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['d*(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))
                                else:
                                    for var in sub_terms_list_copy_1:
                                        if var[-1] == P3_r[-1]:
                                            sub_list_var_index = sub_terms_list_copy_1.index(var)
                                            sub_terms_list_copy_1[sub_list_var_index] = f"{var[:-2]}{right_P3_commutator_index}"
                                            break

                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))

                                # Second term
                                sub_terms_list_copy_2[sub_terms_list_copy_2.index(P3_r)] = f"P3{right_P3_commutator_index}"
                                left_P3_commutator_pos = sub_terms_list_copy_2.index(P3_var_on_left_commute_list)
                                sub_terms_list_copy_2.pop(left_P3_commutator_pos)
                                right_P3_commutator_pos = sub_terms_list_copy_2.index(P3_var)
                                sub_terms_list_copy_2.pop(right_P3_commutator_pos)
                                if P3_var_on_left_commute_list[-1] == P3_r[-1]:
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-d*(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                                else:
                                    for var in sub_terms_list_copy_2:
                                        if var[-1] == P3_r[-1]:
                                            sub_list_var_index = sub_terms_list_copy_2.index(var)
                                            sub_terms_list_copy_2[sub_list_var_index] = f"{var[:-2]}{left_P3_commutator_index}"
                                            break

                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # Now add in the Riemann tensor operator that acts on the fields themselves
                            sub_terms_list_copy_1 = sub_terms_list.copy()
                            sub_terms_list_copy_2 = sub_terms_list.copy()

                            LHS_commutator_P3_index = sub_terms_list_copy_1.index(P3_var_on_left_commute_list)
                            sub_terms_list_copy_1[LHS_commutator_P3_index] = f"a3{left_P3_commutator_index}*U3{right_P3_commutator_index}"
                            RHS_commutator_P3_index = sub_terms_list_copy_1.index(P3_var)
                            sub_terms_list_copy_1.pop(RHS_commutator_P3_index)

                            if sub_terms_list_copy_1[0][0] != '-':
                                sub_terms_list_copy_1 = ['-(1/l^2)'] + sub_terms_list_copy_1
                            else:
                                sub_terms_list_copy_1 = ['(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_1))

                            LHS_commutator_P3_index = sub_terms_list_copy_2.index(P3_var_on_left_commute_list)
                            sub_terms_list_copy_2[LHS_commutator_P3_index] = f"a3{right_P3_commutator_index}*U3{left_P3_commutator_index}"
                            RHS_commutator_P3_index = sub_terms_list_copy_2.index(P3_var)
                            sub_terms_list_copy_2.pop(RHS_commutator_P3_index)

                            if sub_terms_list_copy_2[0][0] != '-':
                                sub_terms_list_copy_2 = ['(1/l^2)'] + sub_terms_list_copy_2
                            else:
                                sub_terms_list_copy_2 = ['-(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # Swap the P3s in both the main list and the P3_to_commute list
                            sub_terms_list[P3_var_on_left_sub_list], sub_terms_list[P3_index_sub_list] = (
                                sub_terms_list[P3_index_sub_list],
                                sub_terms_list[P3_var_on_left_sub_list]
                            )
                            P3_to_commute[P3_index_commute_list - 1], P3_to_commute[P3_index_commute_list] = (
                                P3_to_commute[P3_index_commute_list],
                                P3_to_commute[P3_index_commute_list - 1]
                            )

                        # -------------------------------------------------------
                        # Commute P2s
                        # -------------------------------------------------------
                        P2_to_commute = [var for var in sub_terms_list if var.startswith('P2')]
                        while P2_to_commute.index(P2_var) != 0:
                            P2_index_commute_list = P2_to_commute.index(P2_var)
                            P2_index_sub_list = sub_terms_list.index(P2_var)
                            P2_var_on_left_commute_list = P2_to_commute[P2_index_commute_list - 1]
                            P2_var_on_left_sub_list = sub_terms_list.index(P2_var_on_left_commute_list)
                            left_P2_commutator_index = P2_var_on_left_commute_list[-2:]
                            right_P2_commutator_index = P2_var[-2:]

                            P2s_on_right = [
                                var for var in P2_to_commute
                                if P2_to_commute.index(var) > P2_index_commute_list
                            ]
                            for P2_r in P2s_on_right:
                                sub_terms_list_copy_1 = sub_terms_list.copy()
                                sub_terms_list_copy_2 = sub_terms_list.copy()

                                # First term
                                sub_terms_list_copy_1[sub_terms_list_copy_1.index(P2_r)] = f"P2{left_P2_commutator_index}"
                                left_P2_commutator_pos = sub_terms_list_copy_1.index(P2_var_on_left_commute_list)
                                sub_terms_list_copy_1.pop(left_P2_commutator_pos)
                                right_P2_commutator_pos = sub_terms_list_copy_1.index(P2_var)
                                sub_terms_list_copy_1.pop(right_P2_commutator_pos)
                                if P2_var[-1] == P2_r[-1]:
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['d*(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))

                                else:
                                    for var in sub_terms_list_copy_1:
                                        if var[-1] == P2_r[-1]:
                                            sub_list_var_index = sub_terms_list_copy_1.index(var)
                                            sub_terms_list_copy_1[sub_list_var_index] = f"{var[:-2]}{right_P2_commutator_index}"
                                            break

                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))

                                # Second term
                                sub_terms_list_copy_2[sub_terms_list_copy_2.index(P2_r)] = f"P2{right_P2_commutator_index}"
                                left_P2_commutator_pos = sub_terms_list_copy_2.index(P2_var_on_left_commute_list)
                                sub_terms_list_copy_2.pop(left_P2_commutator_pos)
                                right_P2_commutator_pos = sub_terms_list_copy_2.index(P2_var)
                                sub_terms_list_copy_2.pop(right_P2_commutator_pos)
                                if P2_var_on_left_commute_list[-1] == P2_r[-1]:
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-d*(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                                else:
                                    for var in sub_terms_list_copy_2:
                                        if var[-1] == P2_r[-1]:
                                            sub_list_var_index = sub_terms_list_copy_2.index(var)
                                            sub_terms_list_copy_2[sub_list_var_index] = f"{var[:-2]}{left_P2_commutator_index}"
                                            break

                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # Now we add in the riemann tensor operator that acts on the fields themselves
                            sub_terms_list_copy_1 = sub_terms_list.copy()
                            sub_terms_list_copy_2 = sub_terms_list.copy()

                            LHS_commutator_P2_index = sub_terms_list_copy_1.index(P2_var_on_left_commute_list)
                            sub_terms_list_copy_1[LHS_commutator_P2_index] = f"a2{left_P2_commutator_index}*U2{right_P2_commutator_index}"
                            RHS_commutator_P2_index = sub_terms_list_copy_1.index(P2_var)
                            sub_terms_list_copy_1.pop(RHS_commutator_P2_index)

                            if sub_terms_list_copy_1[0][0] != '-':
                                sub_terms_list_copy_1 = ['-(1/l^2)'] + sub_terms_list_copy_1
                            else:
                                sub_terms_list_copy_1 = ['(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_1))

                            LHS_commutator_P2_index = sub_terms_list_copy_2.index(P2_var_on_left_commute_list)
                            sub_terms_list_copy_2[LHS_commutator_P2_index] = f"a2{right_P2_commutator_index}*U2{left_P2_commutator_index}"
                            RHS_commutator_P2_index = sub_terms_list_copy_2.index(P2_var)
                            sub_terms_list_copy_2.pop(RHS_commutator_P2_index)

                            if sub_terms_list_copy_2[0][0] != '-':
                                sub_terms_list_copy_2 = ['(1/l^2)'] + sub_terms_list_copy_2
                            else:
                                sub_terms_list_copy_2 = ['-(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # Swap the orders of P2s in both lists for the next iteration
                            sub_terms_list[P2_var_on_left_sub_list], sub_terms_list[P2_index_sub_list] = (
                                sub_terms_list[P2_index_sub_list],
                                sub_terms_list[P2_var_on_left_sub_list]
                            )
                            P2_to_commute[P2_index_commute_list - 1], P2_to_commute[P2_index_commute_list] = (
                                P2_to_commute[P2_index_commute_list],
                                P2_to_commute[P2_index_commute_list - 1]
                            )
                        # only want to move one instance - break out of two loops and go to next term
                        break

        if P3P2_present:
            # now that sub_terms_list has the variables in the desired order, we can replace P3P2 with (1/2)*(P1_a*P1^a - P2_a*P2^a - P3_a*P3^a)
            # create copies again to construct desired terms
            sub_terms_list_B1 = sub_terms_list.copy()
            sub_terms_list_B2 = sub_terms_list.copy()
            sub_terms_list_B3 = sub_terms_list.copy()
            # can just use P3_var and P2_var - in first list, delete P3_var and replace the P2_var with the box operator (same index)
            # but need the replaced terms to be on the left of all derivatives due to the partial integration
            P2_var_index = sub_terms_list.index(P2_var)
            leftmost_P_index = P2_var_index
            for var in sub_terms_list_B1:
                # find left most P index
                if var.startswith('P'):
                    leftmost_P_index = sub_terms_list_B1.index(var)
                    break
            # then replace with relevant box operators
            lett_index_to_use = P3_var_stored[-1]
            sub_terms_list_B1.insert(leftmost_P_index, f'P1_{lett_index_to_use}*P1^{lett_index_to_use}')
            sub_terms_list_B2.insert(leftmost_P_index, f'P2_{lett_index_to_use}*P2^{lett_index_to_use}')
            sub_terms_list_B3.insert(leftmost_P_index, f'P3_{lett_index_to_use}*P3^{lett_index_to_use}')
            # B2 term is multiplied by 1/2, while B1 and B3 are multtiplied by -1/2 

            P3_var_index_B1 = sub_terms_list_B1.index(P3_var_stored)
            P3_var_index_B2 = sub_terms_list_B2.index(P3_var_stored)
            P3_var_index_B3 = sub_terms_list_B3.index(P3_var_stored)
            # delete P3_var in all lists
            sub_terms_list_B1.pop(P3_var_index_B1)
            sub_terms_list_B2.pop(P3_var_index_B2)
            sub_terms_list_B3.pop(P3_var_index_B3)

            # and delete the P2_var too
            P2_var_index_B1 = sub_terms_list_B1.index(P2_var)
            P2_var_index_B2 = sub_terms_list_B2.index(P2_var)
            P2_var_index_B3 = sub_terms_list_B3.index(P2_var)
            sub_terms_list_B1.pop(P2_var_index_B1)
            sub_terms_list_B2.pop(P2_var_index_B2)
            sub_terms_list_B3.pop(P2_var_index_B3)
            # B2 term is multiplied by 1/2, while B1 and B3 are multiplied by -1/2 
            if sub_terms_list_B1[0][0] == '-':
                sub_terms_list_B1 = ['-(1/2)'] + [sub_terms_list_B1[0][1:]] + sub_terms_list_B1[1:]
            else:
                sub_terms_list_B1 = ['(1/2)'] + sub_terms_list_B1
            if sub_terms_list_B2[0][0] == '-':
                sub_terms_list_B2 = ['(1/2)'] + [sub_terms_list_B2[0][1:]] + sub_terms_list_B2[1:]
            else:
                sub_terms_list_B2 = ['-(1/2)'] + sub_terms_list_B2
            if sub_terms_list_B3[0][0] == '-':
                sub_terms_list_B3 = ['(1/2)'] + [sub_terms_list_B3[0][1:]] + sub_terms_list_B3[1:]
            else:
                sub_terms_list_B3 = ['-(1/2)'] + sub_terms_list_B3
            
            # then we can construct the terms by joining with *
            expanded_terms.append('*'.join(sub_terms_list_B1))
            expanded_terms.append('*'.join(sub_terms_list_B2))
            expanded_terms.append('*'.join(sub_terms_list_B3))
        else:
            # just reconstruct the original term
            expanded_terms.append('*'.join(sub_terms_list))

    return pull_minus_signs_to_front(' + '.join(expanded_terms))

def pull_P1P2s_left(equation):
    # pulls the leftmost P1P2 to the left to integrate by parts
    # break the equation down into individual terms by top level '+' symbols
    equation = replace_minuses(equation)
    terms = top_level_split(equation)
    expanded_terms = []
    for term in terms:
        # break into individual variables by * signs
        # this will be the main list which will be rearranged as we go
        sub_terms_list = term.split('*')
        # check over sub_terms_list for P1 variables
        P1P2_present = False
        for P1_var in sub_terms_list:
            if P1P2_present == True:
                break
            elif P1_var.startswith('P1'):
                P1_var_stored = P1_var
                # store the letter index and position index of the P1
                P1_lett_index = P1_var[-1]
                P1_pos_index = sub_terms_list.index(P1_var)
                for P2_var in sub_terms_list:
                    if P2_var.startswith('P2') and P2_var[-1] == P1_lett_index:
                        P1P2_present = True
                        # store the P2 index (letter index should be same as P1 if relevant)
                        P2_pos_index = sub_terms_list.index(P2_var)

                        # -------------------------------------------------------
                        # Now we want to move all U2, U1, Z1, Z2, Z3 left of the leftmost P
                        # before we can commute P2s/P1s
                        # -------------------------------------------------------
                        U2_U1_Z1_Z2_Z3_list = [var for var in sub_terms_list if var.startswith('U1') or var.startswith('U2') or var.startswith('U3') or var.startswith('Z2') or var.startswith('Z1') or var.startswith('Z3')]
                        # Create a list with those items removed
                        remaining_list = [
                            var for var in sub_terms_list
                            if var not in U2_U1_Z1_Z2_Z3_list
                        ]

                        # Find the position of the leftmost P in the remaining_list
                        for var in remaining_list:
                            if var.startswith('P'):
                                leftmost_P = var
                                leftmost_P_pos = remaining_list.index(var)
                                break

                        # Insert all U2/U1/Z1/Z2/Z3 vars to the left of that leftmost P
                        for U2_U1_Z1_Z2_Z3 in U2_U1_Z1_Z2_Z3_list:
                            remaining_list.insert(leftmost_P_pos, U2_U1_Z1_Z2_Z3)

                        # Reassign sub_terms_list
                        sub_terms_list = remaining_list

                        # -------------------------------------------------------
                        # Commute P1s
                        # -------------------------------------------------------
                        P1_to_commute = [var for var in sub_terms_list if var.startswith('P1')]
                        while P1_to_commute.index(P1_var) != 0:
                            P1_index_commute_list = P1_to_commute.index(P1_var)
                            P1_index_sub_list = sub_terms_list.index(P1_var)
                            P1_var_on_left_commute_list = P1_to_commute[P1_index_commute_list - 1]
                            P1_var_on_left_sub_list = sub_terms_list.index(P1_var_on_left_commute_list)
                            left_P1_commutator_index = P1_var_on_left_commute_list[-2:]
                            right_P1_commutator_index = P1_var[-2:]

                            # Create terms from Riemann tensor contracting with P1s on the right
                            P1s_on_right = [
                                var for var in P1_to_commute
                                if P1_to_commute.index(var) > P1_index_commute_list
                            ]
                            for P1_r in P1s_on_right:
                                # copy 1
                                sub_terms_list_copy_1 = sub_terms_list.copy()
                                # copy 2
                                sub_terms_list_copy_2 = sub_terms_list.copy()

                                # First term
                                sub_terms_list_copy_1[sub_terms_list_copy_1.index(P1_r)] = f"P1{left_P1_commutator_index}"
                                left_P1_commutator_pos = sub_terms_list_copy_1.index(P1_var_on_left_commute_list)
                                sub_terms_list_copy_1.pop(left_P1_commutator_pos)
                                right_P1_commutator_pos = sub_terms_list_copy_1.index(P1_var)
                                sub_terms_list_copy_1.pop(right_P1_commutator_pos)
                                if P1_var[-1] == P1_r[-1]:
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['d*(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))
                                else:
                                    for var in sub_terms_list_copy_1:
                                        if var[-1] == P1_r[-1]:
                                            sub_list_var_index = sub_terms_list_copy_1.index(var)
                                            sub_terms_list_copy_1[sub_list_var_index] = f"{var[:-2]}{right_P1_commutator_index}"
                                            break

                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))

                                # Second term
                                sub_terms_list_copy_2[sub_terms_list_copy_2.index(P1_r)] = f"P1{right_P1_commutator_index}"
                                left_P1_commutator_pos = sub_terms_list_copy_2.index(P1_var_on_left_commute_list)
                                sub_terms_list_copy_2.pop(left_P1_commutator_pos)
                                right_P1_commutator_pos = sub_terms_list_copy_2.index(P1_var)
                                sub_terms_list_copy_2.pop(right_P1_commutator_pos)
                                if P1_var_on_left_commute_list[-1] == P1_r[-1]:
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-d*(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                                else:
                                    for var in sub_terms_list_copy_2:
                                        if var[-1] == P1_r[-1]:
                                            sub_list_var_index = sub_terms_list_copy_2.index(var)
                                            sub_terms_list_copy_2[sub_list_var_index] = f"{var[:-2]}{left_P1_commutator_index}"
                                            break

                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # Now add in the Riemann tensor operator that acts on the fields themselves
                            sub_terms_list_copy_1 = sub_terms_list.copy()
                            sub_terms_list_copy_2 = sub_terms_list.copy()

                            LHS_commutator_P1_index = sub_terms_list_copy_1.index(P1_var_on_left_commute_list)
                            sub_terms_list_copy_1[LHS_commutator_P1_index] = f"a1{left_P1_commutator_index}*U1{right_P1_commutator_index}"
                            RHS_commutator_P1_index = sub_terms_list_copy_1.index(P1_var)
                            sub_terms_list_copy_1.pop(RHS_commutator_P1_index)

                            if sub_terms_list_copy_1[0][0] != '-':
                                sub_terms_list_copy_1 = ['-(1/l^2)'] + sub_terms_list_copy_1
                            else:
                                sub_terms_list_copy_1 = ['(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_1))

                            LHS_commutator_P1_index = sub_terms_list_copy_2.index(P1_var_on_left_commute_list)
                            sub_terms_list_copy_2[LHS_commutator_P1_index] = f"a1{right_P1_commutator_index}*U1{left_P1_commutator_index}"
                            RHS_commutator_P1_index = sub_terms_list_copy_2.index(P1_var)
                            sub_terms_list_copy_2.pop(RHS_commutator_P1_index)

                            if sub_terms_list_copy_2[0][0] != '-':
                                sub_terms_list_copy_2 = ['(1/l^2)'] + sub_terms_list_copy_2
                            else:
                                sub_terms_list_copy_2 = ['-(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # Swap the P1s in both the main list and the P1_to_commute list
                            sub_terms_list[P1_var_on_left_sub_list], sub_terms_list[P1_index_sub_list] = (
                                sub_terms_list[P1_index_sub_list],
                                sub_terms_list[P1_var_on_left_sub_list]
                            )
                            P1_to_commute[P1_index_commute_list - 1], P1_to_commute[P1_index_commute_list] = (
                                P1_to_commute[P1_index_commute_list],
                                P1_to_commute[P1_index_commute_list - 1]
                            )

                        # -------------------------------------------------------
                        # Commute P2s
                        # -------------------------------------------------------
                        P2_to_commute = [var for var in sub_terms_list if var.startswith('P2')]
                        while P2_to_commute.index(P2_var) != 0:
                            P2_index_commute_list = P2_to_commute.index(P2_var)
                            P2_index_sub_list = sub_terms_list.index(P2_var)
                            P2_var_on_left_commute_list = P2_to_commute[P2_index_commute_list - 1]
                            P2_var_on_left_sub_list = sub_terms_list.index(P2_var_on_left_commute_list)
                            left_P2_commutator_index = P2_var_on_left_commute_list[-2:]
                            right_P2_commutator_index = P2_var[-2:]

                            P2s_on_right = [
                                var for var in P2_to_commute
                                if P2_to_commute.index(var) > P2_index_commute_list
                            ]
                            for P2_r in P2s_on_right:
                                sub_terms_list_copy_1 = sub_terms_list.copy()
                                sub_terms_list_copy_2 = sub_terms_list.copy()

                                # First term
                                sub_terms_list_copy_1[sub_terms_list_copy_1.index(P2_r)] = f"P2{left_P2_commutator_index}"
                                left_P2_commutator_pos = sub_terms_list_copy_1.index(P2_var_on_left_commute_list)
                                sub_terms_list_copy_1.pop(left_P2_commutator_pos)
                                right_P2_commutator_pos = sub_terms_list_copy_1.index(P2_var)
                                sub_terms_list_copy_1.pop(right_P2_commutator_pos)
                                if P2_var[-1] == P2_r[-1]:
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['d*(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))
                                else:
                                    for var in sub_terms_list_copy_1:
                                        if var[-1] == P2_r[-1]:
                                            sub_list_var_index = sub_terms_list_copy_1.index(var)
                                            sub_terms_list_copy_1[sub_list_var_index] = f"{var[:-2]}{right_P2_commutator_index}"
                                            break

                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))

                                # Second term
                                sub_terms_list_copy_2[sub_terms_list_copy_2.index(P2_r)] = f"P2{right_P2_commutator_index}"
                                left_P2_commutator_pos = sub_terms_list_copy_2.index(P2_var_on_left_commute_list)
                                sub_terms_list_copy_2.pop(left_P2_commutator_pos)
                                right_P2_commutator_pos = sub_terms_list_copy_2.index(P2_var)
                                sub_terms_list_copy_2.pop(right_P2_commutator_pos)
                                if P2_var_on_left_commute_list[-1] == P2_r[-1]:
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-d*(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                                else:
                                    for var in sub_terms_list_copy_2:
                                        if var[-1] == P2_r[-1]:
                                            sub_list_var_index = sub_terms_list_copy_2.index(var)
                                            sub_terms_list_copy_2[sub_list_var_index] = f"{var[:-2]}{left_P2_commutator_index}"
                                            break

                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # Now we add in the riemann tensor operator that acts on the fields themselves
                            sub_terms_list_copy_1 = sub_terms_list.copy()
                            sub_terms_list_copy_2 = sub_terms_list.copy()

                            LHS_commutator_P2_index = sub_terms_list_copy_1.index(P2_var_on_left_commute_list)
                            sub_terms_list_copy_1[LHS_commutator_P2_index] = f"a2{left_P2_commutator_index}*U2{right_P2_commutator_index}"
                            RHS_commutator_P2_index = sub_terms_list_copy_1.index(P2_var)
                            sub_terms_list_copy_1.pop(RHS_commutator_P2_index)

                            if sub_terms_list_copy_1[0][0] != '-':
                                sub_terms_list_copy_1 = ['-(1/l^2)'] + sub_terms_list_copy_1
                            else:
                                sub_terms_list_copy_1 = ['(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_1))

                            LHS_commutator_P2_index = sub_terms_list_copy_2.index(P2_var_on_left_commute_list)
                            sub_terms_list_copy_2[LHS_commutator_P2_index] = f"a2{right_P2_commutator_index}*U2{left_P2_commutator_index}"
                            RHS_commutator_P2_index = sub_terms_list_copy_2.index(P2_var)
                            sub_terms_list_copy_2.pop(RHS_commutator_P2_index)

                            if sub_terms_list_copy_2[0][0] != '-':
                                sub_terms_list_copy_2 = ['(1/l^2)'] + sub_terms_list_copy_2
                            else:
                                sub_terms_list_copy_2 = ['-(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # Swap the orders of P2s in both lists for the next iteration
                            sub_terms_list[P2_var_on_left_sub_list], sub_terms_list[P2_index_sub_list] = (
                                sub_terms_list[P2_index_sub_list],
                                sub_terms_list[P2_var_on_left_sub_list]
                            )
                            P2_to_commute[P2_index_commute_list - 1], P2_to_commute[P2_index_commute_list] = (
                                P2_to_commute[P2_index_commute_list],
                                P2_to_commute[P2_index_commute_list - 1]
                            )
                        # only want to move one instance - break out of two loops and go to next term
                        break

        if P1P2_present:
            # now that sub_terms_list has the variables in the desired order,
            # we can replace P1P2 with (1/2)*(P3_a*P3^a - P2_a*P2^a - P1_a*P1^a) 
            sub_terms_list_B1 = sub_terms_list.copy()
            sub_terms_list_B2 = sub_terms_list.copy()
            sub_terms_list_B3 = sub_terms_list.copy()
            # can just use P1_var and P2_var - in first list, delete the P1_var and replace the P2_var with the box operator (same index)
            # but need the replaced terms to be on the left of all derivatives due to the partial integration
            P2_var_index = sub_terms_list.index(P2_var)
            leftmost_P_index = P2_var_index
            for var in sub_terms_list_B1:
                # find left most P index
                if var.startswith('P'):
                    leftmost_P_index = sub_terms_list_B1.index(var)
                    break
            # then replace with relevant box operators
            lett_index_to_use = P2_var[-1]
            sub_terms_list_B1.insert(leftmost_P_index, f'P1_{lett_index_to_use}*P1^{lett_index_to_use}')
            sub_terms_list_B2.insert(leftmost_P_index, f'P2_{lett_index_to_use}*P2^{lett_index_to_use}')
            sub_terms_list_B3.insert(leftmost_P_index, f'P3_{lett_index_to_use}*P3^{lett_index_to_use}')
            # B2 term is multiplied by 1/2, while B1 and B3 are multtiplied by -1/2 

            P1_var_index_B1 = sub_terms_list_B1.index(P1_var_stored)
            P1_var_index_B2 = sub_terms_list_B2.index(P1_var_stored)
            P1_var_index_B3 = sub_terms_list_B3.index(P1_var_stored)
            # delete P1_var in all lists
            sub_terms_list_B1.pop(P1_var_index_B1)
            sub_terms_list_B2.pop(P1_var_index_B2)
            sub_terms_list_B3.pop(P1_var_index_B3)

            # and delete the P2_var too
            P2_var_index_B1 = sub_terms_list_B1.index(P2_var)
            P2_var_index_B2 = sub_terms_list_B2.index(P2_var)
            P2_var_index_B3 = sub_terms_list_B3.index(P2_var)
            sub_terms_list_B1.pop(P2_var_index_B1)
            sub_terms_list_B2.pop(P2_var_index_B2)
            sub_terms_list_B3.pop(P2_var_index_B3)
            # B2 term is multiplied by 1/2, while B1 and B3 are multiplied by -1/2 
            if sub_terms_list_B1[0][0] == '-':
                sub_terms_list_B1 = ['(1/2)'] + [sub_terms_list_B1[0][1:]] + sub_terms_list_B1[1:]
            else:
                sub_terms_list_B1 = ['-(1/2)'] + sub_terms_list_B1
            if sub_terms_list_B2[0][0] == '-':
                sub_terms_list_B2 = ['(1/2)'] + [sub_terms_list_B2[0][1:]] + sub_terms_list_B2[1:]
            else:
                sub_terms_list_B2 = ['-(1/2)'] + sub_terms_list_B2
            if sub_terms_list_B3[0][0] == '-':
                sub_terms_list_B3 = ['-(1/2)'] + [sub_terms_list_B3[0][1:]] + sub_terms_list_B3[1:]
            else:
                sub_terms_list_B3 = ['(1/2)'] + sub_terms_list_B3
            
            # then we can construct the terms by joining with '*'
            expanded_terms.append('*'.join(sub_terms_list_B1))
            expanded_terms.append('*'.join(sub_terms_list_B2))
            expanded_terms.append('*'.join(sub_terms_list_B3))
        else:
            # just reconstruct the original term
            expanded_terms.append('*'.join(sub_terms_list))

    return pull_minus_signs_to_front(' + '.join(expanded_terms))

def pull_all_PiPjs_left(equation):
    return move_all_a_left(pull_P1P2s_left(move_all_a_left(pull_P3P1s_left(move_all_a_left(pull_P3P2s_left(equation))))))

def pull_B1s_right(equation):
    # pulls the leftmost B1 to the right to impose on shell conditions
    # break the equation down into individual terms by top level '+' symbols
    equation = replace_minuses(equation)
    terms = top_level_split(equation)
    expanded_terms = []
    for term in terms:
        # break into individual variables by * signs
        # this will be the main list which will be rearranged as we go
        sub_terms_list = term.split('*')
        # check over sub_terms_list for P1 variables
        B1_present = False
        for P1_var in reversed(sub_terms_list):
            if B1_present == True:
                break
            elif P1_var.startswith('P1'):
                P1_var_stored = P1_var
                # store the letter index and position index of the P1
                P1_lett_index = P1_var[-1]
                P1_pos_index = sub_terms_list.index(P1_var)
                for P11_var in reversed(sub_terms_list):
                    if P11_var.startswith('P1') and P11_var[-1] == P1_lett_index and P11_var != P1_var:
                        B1_present = True
                        # store the P11 index (letter index should be same as P1 if relevant)
                        P11_pos_index = sub_terms_list.index(P11_var)

                        # -------------------------------------------------------
                        # Now we want to move all U1, Z2, Z3 left of the leftmost P
                        # before we can commute P1s
                        # -------------------------------------------------------
                        U1_Z2_Z3_list = [var for var in sub_terms_list if var.startswith('U1') or var.startswith('U2') or var.startswith('U3') or var.startswith('Z2') or var.startswith('Z1') or var.startswith('Z3')]
                        # Create a list with those items removed
                        remaining_list = [
                            var for var in sub_terms_list
                            if var not in U1_Z2_Z3_list
                        ]

                        # Find the position of the leftmost P in the remaining_list
                        for var in remaining_list:
                            if var.startswith('P'):
                                leftmost_P = var
                                leftmost_P_pos = remaining_list.index(var)
                                break

                        # Insert all U1/Z2/Z3 vars to the left of that leftmost P
                        for U1_Z2_Z3 in U1_Z2_Z3_list:
                            remaining_list.insert(leftmost_P_pos, U1_Z2_Z3)

                        # Reassign sub_terms_list
                        sub_terms_list = remaining_list

                        # -------------------------------------------------------
                        # Commute P1s
                        # -------------------------------------------------------
                        P1_to_commute = [var for var in sub_terms_list if var.startswith('P1')]
                        while P1_to_commute.index(P1_var) != (len(P1_to_commute)-1):
                            # find P1 index in the to_commute list
                            P1_index_commute_list = P1_to_commute.index(P1_var)
                            # find P1 index in the sub_terms_list
                            P1_index_sub_list = sub_terms_list.index(P1_var)
                            # store the P1 variable on the next right to be commuted through
                            P1_var_on_right_commute_list = P1_to_commute[P1_index_commute_list+1]
                            # find this variable's position in the sub_terms_list
                            P1_var_on_right_sub_list = sub_terms_list.index(P1_var_on_right_commute_list)
                            # extract relevant information about the P3s being commuted (indices)
                            # used for constructing the indices in the riemann tensor
                            right_P1_commutator_index = P1_var_on_right_commute_list[-2:]
                            left_P1_commutator_index = P1_var[-2:]
                            # create terms arising from the riemann tensor contracting with P1s on the right
                            P1s_on_right = [var for var in P1_to_commute if P1_to_commute.index(var) > P1_index_commute_list+1]
                            # if letter index on P1 being commuted through is the same, can just swap the order without concern
                            if P1_var_on_right_commute_list[-1] == P1_var[-1]:
                                sub_terms_list[P1_var_on_right_sub_list], sub_terms_list[P1_index_sub_list] = sub_terms_list[P1_index_sub_list], sub_terms_list[P1_var_on_right_sub_list]
                                P1_to_commute[P1_index_commute_list+1], P1_to_commute[P1_index_commute_list] = P1_to_commute[P1_index_commute_list], P1_to_commute[P1_index_commute_list+1]
                                continue
                            for P1_r in P1s_on_right:
                                # create copy to manipulate terms without altering master list
                                sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*P1_a*g_cb
                                sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*P1_b*g_ca
                                # construct the first term using the first copy
                                # replace P1_r in the copy list with P1_{index from left P1 in commutator}
                                sub_terms_list_copy_1[sub_terms_list_copy_1.index(P1_r)] = f"P1{left_P1_commutator_index}"
                                # then delete the P1 that was in the LHS of the commutator the copy list since this is accounted for in the commutator
                                left_P1_commutator_pos = sub_terms_list_copy_1.index(P1_var_on_right_commute_list)
                                sub_terms_list_copy_1.pop(left_P1_commutator_pos)
                                # then also delete the P1 that was in the RHS of the commutator from the list
                                right_P1_commutator_pos = sub_terms_list_copy_1.index(P1_var)
                                sub_terms_list_copy_1.pop(right_P1_commutator_pos)
                                # then in this term, we will also have a metric to contract: g_{P1_r index}_{index from right P1 in commutator}
                                # also could have P1_r index = index from right P1 in commutator -> just include a factor of d (dimension)
                                if P1_var_on_right_commute_list[-1] == P1_r[-1]:
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['d*(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))
                                else:
                                    # so find another term with one of these indices (choose P1_r cos why not) and change its letter index accordingly
                                    for var in sub_terms_list_copy_1:
                                        # check if index is same as P1_r -> this var will be contracted with the metric
                                        if var[-1] == P1_r[-1]:
                                            # find the index of this var in the list copy
                                            sub_list_var_index = sub_terms_list_copy_1.index(var)
                                            # if the P1_r index is being contracted with var, then the remaining index in the index of the
                                            # P1 in the RHS of the commutator
                                            sub_terms_list_copy_1[sub_list_var_index] = f"{var[:-2]}{right_P1_commutator_index}"
                                            break
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))

                                # now repeat the process but for the second term; same structure, except for a sign and swap of indices
                                # construct the second term using the second copy
                                # delete the P1_r in the copy list since this is accounted for in the commutator
                                left_P1_commutator_pos = sub_terms_list_copy_2.index(P1_var_on_right_commute_list)
                                sub_terms_list_copy_2.pop(left_P1_commutator_pos)
                                # then also delete the P1 that was in the RHS of the commutator from the list
                                right_P3_commutator_pos = sub_terms_list_copy_2.index(P1_var)
                                sub_terms_list_copy_2.pop(right_P1_commutator_pos)
                                # replace P1_r in the copy list with P1_{index from right P1 in commutator}
                                sub_terms_list_copy_2[sub_terms_list_copy_2.index(P1_r)] = f"P1{right_P1_commutator_index}"
                                # then in this term, we will also have a metric to contract: g_{P1_r index}_{index from left P1 in commutator}
                                # also could have P1_r index = index from right P1 in commutator -> just include a factor of d (dimension)
                                if P1_var[-1] == P1_r[-1]:
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['-d*(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))
                                else:
                                    # so find another term with one of these indices (choose P1_r cos why not) and change its letter index accordingly
                                    for var in sub_terms_list_copy_2:
                                        # check if index is same as P1_r -> this var will be contracted with the metric
                                        if var[-1] == P1_r[-1]:
                                            # find the index of this var in the list copy
                                            sub_list_var_index = sub_terms_list_copy_2.index(var)
                                            # if the P3_r index is being contracted with var, then the remaining index in the index of the
                                            # P3 in the LHS of the commutator
                                            sub_terms_list_copy_2[sub_list_var_index] = f"{var[:-2]}{left_P1_commutator_index}"
                                            break
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                            
                            # now we have to add in the riemann tensor operator that acts on the fields themselves.
                            # again, make two separate terms
                            sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*a1_b*U1_a
                            sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*a1_a*U1_b
                            # the indices are now only taken from P_var and P_var_on_left
                            # replace the RHS commutator P1 term with aU term
                            LHS_commutator_P1_index = sub_terms_list_copy_1.index(P1_var_on_right_commute_list)
                            sub_terms_list_copy_1[LHS_commutator_P1_index] = f"a1{left_P1_commutator_index}*U1{right_P1_commutator_index}"
                            # then delete the RHS commutator P1 term from the list
                            RHS_commutator_P1_index = sub_terms_list_copy_1.index(P1_var)
                            sub_terms_list_copy_1.pop(RHS_commutator_P1_index)
                            # then add in the 1/l^2 term, accounting for signs
                            if sub_terms_list_copy_1[0][0] != '-':
                                sub_terms_list_copy_1 = ['-(1/l^2)'] + sub_terms_list_copy_1
                            else:
                                sub_terms_list_copy_1 = ['(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_1))

                            # then the second term
                            # replace the LHS commutator P1 term with aU term
                            LHS_commutator_P1_index = sub_terms_list_copy_2.index(P1_var_on_right_commute_list)
                            sub_terms_list_copy_2[LHS_commutator_P1_index] = f"a1{right_P1_commutator_index}*U1{left_P1_commutator_index}"
                            # then delete the RHS commutator P1 term from the list
                            RHS_commutator_P1_index = sub_terms_list_copy_2.index(P1_var)
                            sub_terms_list_copy_2.pop(RHS_commutator_P1_index)
                            # then add in the 1/l^2 term, accounting for signs
                            if sub_terms_list_copy_2[0][0] != '-':
                                sub_terms_list_copy_2 = ['(1/l^2)'] + sub_terms_list_copy_2
                            else:
                                sub_terms_list_copy_2 = ['-(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_2))
                            # lastly swap the orders of P1s in both lists for the next iteration
                            sub_terms_list[P1_var_on_right_sub_list], sub_terms_list[P1_index_sub_list] = sub_terms_list[P1_index_sub_list], sub_terms_list[P1_var_on_right_sub_list]
                            P1_to_commute[P1_index_commute_list+1], P1_to_commute[P1_index_commute_list] = P1_to_commute[P1_index_commute_list], P1_to_commute[P1_index_commute_list+1]
                        # -------------------------------------------------------
                        # Commute second P1
                        # -------------------------------------------------------
                        P11_to_commute = [var for var in sub_terms_list if var.startswith('P1')]
                        while P11_to_commute.index(P11_var) != (len(P11_to_commute)-1):
                            # find P11 index in the to_commute list
                            P11_index_commute_list = P11_to_commute.index(P11_var)
                            # find P11 index in the sub_terms_list
                            P11_index_sub_list = sub_terms_list.index(P11_var)
                            # store the P11 variable on the next right to be commuted through
                            P11_var_on_right_commute_list = P11_to_commute[P11_index_commute_list+1]
                            # find this variable's position in the sub_terms_list
                            P11_var_on_right_sub_list = sub_terms_list.index(P11_var_on_right_commute_list)
                            # extract relevant information about the P1s being commuted (indices)
                            # used for constructing the indices in the riemann tensor
                            right_P11_commutator_index = P11_var_on_right_commute_list[-2:]
                            left_P11_commutator_index = P11_var[-2:]
                            # create terms arising from the riemann tensor contracting with P1s on the right
                            P11s_on_right = [var for var in P11_to_commute if P11_to_commute.index(var) > P11_index_commute_list+1]
                            # if letter index on P1 being commuted through is the same, can just swap the order without concern
                            if P11_var_on_right_commute_list[-1] == P11_var[-1]:
                                sub_terms_list[P11_var_on_right_sub_list], sub_terms_list[P11_index_sub_list] = sub_terms_list[P11_index_sub_list], sub_terms_list[P11_var_on_right_sub_list]
                                P11_to_commute[P11_index_commute_list+1], P11_to_commute[P11_index_commute_list] = P11_to_commute[P11_index_commute_list], P11_to_commute[P11_index_commute_list+1]
                                continue
                            for P11_r in P11s_on_right:
                                # create copy to manipulate terms without altering master list
                                sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*P1_a*g_cb
                                sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*P1_b*g_ca
                                # construct the first term using the first copy
                                # then delete the P1 that was in the LHS of the commutator from the copy list since this is accounted 
                                # for in the commutator
                                left_P11_commutator_pos = sub_terms_list_copy_1.index(P11_var_on_right_commute_list)
                                sub_terms_list_copy_1.pop(left_P11_commutator_pos)
                                # then also delete the P1 that was in the RHS of the commutator from the list
                                right_P11_commutator_pos = sub_terms_list_copy_1.index(P11_var)
                                sub_terms_list_copy_1.pop(right_P11_commutator_pos)
                                # replace P1_r in the copy list with P1_{index from left P1 in commutator}
                                sub_terms_list_copy_1[sub_terms_list_copy_1.index(P11_r)] = f"P1{left_P11_commutator_index}"
                                # then in this term, we will also have a metric to contract: g_{P11_r index}_{index from right P1 in commutator}
                                # also could have P1_r index = index from right P1 in commutator -> just include a factor of d (dimension)
                                if P11_var_on_right_commute_list[-1] == P11_r[-1]:
                                    # now add in the d*1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['d*(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))
                                else:
                                    # so find another term with one of these indices (choose P1_r cos why not) and change its letter index accordingly
                                    for var in sub_terms_list_copy_1:
                                        # check if index is same as P11_r -> this var will be contracted with the metric
                                        if var[-1] == P11_r[-1]:
                                            # find the index of this var in the list copy
                                            sub_list_var_index = sub_terms_list_copy_1.index(var)
                                            # if the P1_r index is being contracted with var, then the remaining index in the index of the
                                            # P1 in the RHS of the commutator
                                            sub_terms_list_copy_1[sub_list_var_index] = f"{var[:-2]}{right_P11_commutator_index}"
                                            break
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))

                                # now repeat the process but for the second term; same structure, except for a sign and swap of indices
                                # construct the second term using the second copy
                                # delete the P11_r in the copy list since this is accounted for in the commutator
                                left_P11_commutator_pos = sub_terms_list_copy_2.index(P11_var_on_right_commute_list)
                                sub_terms_list_copy_2.pop(left_P11_commutator_pos)
                                # then also delete the P1 that was in the RHS of the commutator from the list
                                right_P3_commutator_pos = sub_terms_list_copy_2.index(P11_var)
                                sub_terms_list_copy_2.pop(right_P11_commutator_pos)
                                # replace P11_r in the copy list with P1_{index from right P1 in commutator}
                                sub_terms_list_copy_2[sub_terms_list_copy_2.index(P11_r)] = f"P1{right_P11_commutator_index}"
                                # then in this term, we will also have a metric to contract: g_{P1_r index}_{index from left P1 in commutator}
                                # also could have P1_r index = index from right P1 in commutator -> just include a factor of d (dimension)
                                if P11_var[-1] == P11_r[-1]:
                                    # now add in the d*1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-d*(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['d*(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                                else:
                                    # so find another term with one of these indices (choose P11_r cos why not) and change its letter index accordingly
                                    for var in sub_terms_list_copy_2:
                                        # check if index is same as P11_r -> this var will be contracted with the metric
                                        if var[-1] == P11_r[-1]:
                                            # find the index of this var in the list copy
                                            sub_list_var_index = sub_terms_list_copy_2.index(var)
                                            # if the P11_r index is being contracted with var, then the remaining index in the index of the
                                            # P1 in the LHS of the commutator
                                            sub_terms_list_copy_2[sub_list_var_index] = f"{var[:-2]}{left_P11_commutator_index}"
                                            break
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                            
                            # now we have to add in the riemann tensor operator that acts on the fields themselves.
                            # again, make two separate terms
                            sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*a1_b*U1_a
                            sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*a1_a*U1_b
                            # the indices are now only taken from P11_var and P11_var_on_right
                            # replace the LHS commutator P1 term with aU term
                            LHS_commutator_P11_index = sub_terms_list_copy_1.index(P11_var_on_right_commute_list)
                            sub_terms_list_copy_1[LHS_commutator_P11_index] = f"a1{left_P11_commutator_index}*U1{right_P11_commutator_index}"
                            # then delete the RHS commutator P1 term from the list
                            RHS_commutator_P11_index = sub_terms_list_copy_1.index(P11_var)
                            sub_terms_list_copy_1.pop(RHS_commutator_P11_index)
                            # then add in the 1/l^2 term, accounting for signs
                            if sub_terms_list_copy_1[0][0] != '-':
                                sub_terms_list_copy_1 = ['-(1/l^2)'] + sub_terms_list_copy_1
                            else:
                                sub_terms_list_copy_1 = ['(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_1))

                            # then the second term
                            # replace the LHS commutator P1 term with aU term
                            LHS_commutator_P11_index = sub_terms_list_copy_2.index(P11_var_on_right_commute_list)
                            sub_terms_list_copy_2[LHS_commutator_P11_index] = f"a1{right_P11_commutator_index}*U1{left_P11_commutator_index}"
                            # then delete the RHS commutator P1 term from the list
                            RHS_commutator_P11_index = sub_terms_list_copy_2.index(P11_var)
                            sub_terms_list_copy_2.pop(RHS_commutator_P11_index)
                            # then add in the 1/l^2 term, accounting for signs
                            if sub_terms_list_copy_2[0][0] != '-':
                                sub_terms_list_copy_2 = ['(1/l^2)'] + sub_terms_list_copy_2
                            else:
                                sub_terms_list_copy_2 = ['-(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # lastly swap the orders of P1s in both lists for the next iteration
                            sub_terms_list[P11_var_on_right_sub_list], sub_terms_list[P11_index_sub_list] = sub_terms_list[P11_index_sub_list], sub_terms_list[P11_var_on_right_sub_list]
                            P11_to_commute[P11_index_commute_list+1], P11_to_commute[P11_index_commute_list] = P11_to_commute[P11_index_commute_list], P11_to_commute[P11_index_commute_list+1]
                        # only want to move one instance - break out of two loops and go to next term
                        break

        if B1_present:
            # now that sub_terms_list has the variables in the desired order,
            # we can replace B1 with m1 
            sub_terms_list_B1 = sub_terms_list.copy()
            # can just use P1_var and P11_var - in first list, delete the P1_var and replace the P11_var with the box operator (same index)
            P1_var_index = sub_terms_list.index(P1_var_stored)
            # delete P1_var in list
            sub_terms_list_B1.pop(P1_var_index)
            # find P11_var index
            P11_var_index_B1 = sub_terms_list_B1.index(P11_var)
            # and also delete
            sub_terms_list_B1.pop(P11_var_index_B1)
            # have to include extra on shell constant m1
            if sub_terms_list_B1[0][0] == '-':
                sub_terms_list_B1 = ['-m1'] + [sub_terms_list_B1[0][1:]] + sub_terms_list_B1[1:]
            else:
                sub_terms_list_B1 = ['m1'] + sub_terms_list_B1
            
            # then we can construct the terms by joining with '*'
            expanded_terms.append('*'.join(sub_terms_list_B1))
        else:
            # just reconstruct the original term
            expanded_terms.append('*'.join(sub_terms_list))

    return pull_minus_signs_to_front(' + '.join(expanded_terms))

def pull_B2s_right(equation):
    # pulls the leftmost B2 to the right to impose on shell conditions
    # break the equation down into individual terms by top level '+' symbols
    equation = replace_minuses(equation)
    terms = top_level_split(equation)
    expanded_terms = []
    for term in terms:
        # break into individual variables by * signs
        # this will be the main list which will be rearranged as we go
        sub_terms_list = term.split('*')
        # check over sub_terms_list for P2 variables
        B2_present = False
        for P2_var in reversed(sub_terms_list):
            if B2_present == True:
                break
            elif P2_var.startswith('P2'):
                P2_var_stored = P2_var
                # store the letter index and position index of the P2
                P2_lett_index = P2_var[-1]
                P2_pos_index = sub_terms_list.index(P2_var)
                for P22_var in reversed(sub_terms_list):
                    if P22_var.startswith('P2') and P22_var[-1] == P2_lett_index and P22_var != P2_var:
                        B2_present = True
                        # store the P22 index (letter index should be same as P2 if relevant)
                        P22_pos_index = sub_terms_list.index(P22_var)

                        # -------------------------------------------------------
                        # Now we want to move all U2, Z1, Z3 left of the leftmost P
                        # before we can commute P2s
                        # -------------------------------------------------------
                        U2_Z1_Z3_list = [var for var in sub_terms_list if var.startswith('U1') or var.startswith('U2') or var.startswith('U3') or var.startswith('Z2') or var.startswith('Z1') or var.startswith('Z3')]
                        # Create a list with those items removed
                        remaining_list = [
                            var for var in sub_terms_list
                            if var not in U2_Z1_Z3_list
                        ]

                        # Find the position of the leftmost P in the remaining_list
                        for var in remaining_list:
                            if var.startswith('P'):
                                leftmost_P = var
                                leftmost_P_pos = remaining_list.index(var)
                                break

                        # Insert all U2/Z1/Z3 vars to the left of that leftmost P
                        for U2_Z1_Z3 in U2_Z1_Z3_list:
                            remaining_list.insert(leftmost_P_pos, U2_Z1_Z3)

                        # Reassign sub_terms_list
                        sub_terms_list = remaining_list

                        # -------------------------------------------------------
                        # Commute P2s
                        # -------------------------------------------------------
                        P2_to_commute = [var for var in sub_terms_list if var.startswith('P2')]
                        while P2_to_commute.index(P2_var) != (len(P2_to_commute) - 1):
                            # find P2 index in the to_commute list
                            P2_index_commute_list = P2_to_commute.index(P2_var)
                            # find P2 index in the sub_terms_list
                            P2_index_sub_list = sub_terms_list.index(P2_var)
                            # store the P2 variable on the next right to be commuted through
                            P2_var_on_right_commute_list = P2_to_commute[P2_index_commute_list + 1]
                            # find this variable's position in the sub_terms_list
                            P2_var_on_right_sub_list = sub_terms_list.index(P2_var_on_right_commute_list)
                            # extract relevant information about the P3s being commuted (indices)
                            # used for constructing the indices in the riemann tensor
                            right_P2_commutator_index = P2_var_on_right_commute_list[-2:]
                            left_P2_commutator_index = P2_var[-2:]
                            # create terms arising from the riemann tensor contracting with P2s on the right
                            P2s_on_right = [
                                var for var in P2_to_commute
                                if P2_to_commute.index(var) > P2_index_commute_list + 1
                            ]
                            # if letter index on P2 being commuted through is the same, can just swap the order without concern
                            if P2_var_on_right_commute_list[-1] == P2_var[-1]:
                                sub_terms_list[P2_var_on_right_sub_list], sub_terms_list[P2_index_sub_list] = (
                                    sub_terms_list[P2_index_sub_list],
                                    sub_terms_list[P2_var_on_right_sub_list]
                                )
                                P2_to_commute[P2_index_commute_list + 1], P2_to_commute[P2_index_commute_list] = (
                                    P2_to_commute[P2_index_commute_list],
                                    P2_to_commute[P2_index_commute_list + 1]
                                )
                                continue
                            for P2_r in P2s_on_right:
                                # create copy to manipulate terms without altering master list
                                sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*P2_a*g_cb
                                sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*P2_b*g_ca
                                # construct the first term using the first copy
                                # replace P2_r in the copy list with P2_{index from left P2 in commutator}
                                sub_terms_list_copy_1[sub_terms_list_copy_1.index(P2_r)] = f"P2{left_P2_commutator_index}"
                                # then delete the P2 that was in the LHS of the commutator the copy list since this is accounted for in the commutator
                                left_P2_commutator_pos = sub_terms_list_copy_1.index(P2_var_on_right_commute_list)
                                sub_terms_list_copy_1.pop(left_P2_commutator_pos)
                                # then also delete the P2 that was in the RHS of the commutator from the list
                                right_P2_commutator_pos = sub_terms_list_copy_1.index(P2_var)
                                sub_terms_list_copy_1.pop(right_P2_commutator_pos)
                                # then in this term, we will also have a metric to contract: g_{P2_r index}_{index from right P2 in commutator}
                                # also could have P2_r index = index from right P2 in commutator -> just include a factor of d (dimension)
                                if P2_var_on_right_commute_list[-1] == P2_r[-1]:
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['d*(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))
                                else:
                                    # so find another term with one of these indices (choose P2_r cos why not) and change its letter index accordingly
                                    for var in sub_terms_list_copy_1:
                                        # check if index is same as P2_r -> this var will be contracted with the metric
                                        if var[-1] == P2_r[-1]:
                                            sub_list_var_index = sub_terms_list_copy_1.index(var)
                                            sub_terms_list_copy_1[sub_list_var_index] = f"{var[:-2]}{right_P2_commutator_index}"
                                            break
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))

                                # now repeat the process but for the second term; same structure, except for a sign and swap of indices
                                # construct the second term using the second copy
                                left_P2_commutator_pos = sub_terms_list_copy_2.index(P2_var_on_right_commute_list)
                                sub_terms_list_copy_2.pop(left_P2_commutator_pos)
                                right_P3_commutator_pos = sub_terms_list_copy_2.index(P2_var)
                                sub_terms_list_copy_2.pop(right_P3_commutator_pos)
                                sub_terms_list_copy_2[sub_terms_list_copy_2.index(P2_r)] = f"P2{right_P2_commutator_index}"
                                if P2_var[-1] == P2_r[-1]:
                                    # now add in the 1/l^2 term, accounting for possible minus signs
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['-d*(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))
                                else:
                                    for var in sub_terms_list_copy_2:
                                        if var[-1] == P2_r[-1]:
                                            sub_list_var_index = sub_terms_list_copy_2.index(var)
                                            sub_terms_list_copy_2[sub_list_var_index] = f"{var[:-2]}{left_P2_commutator_index}"
                                            break
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                            
                            # now we have to add in the riemann tensor operator that acts on the fields themselves.
                            # again, make two separate terms
                            sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*a1_b*U2_a
                            sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*a1_a*U2_b
                            LHS_commutator_P2_index = sub_terms_list_copy_1.index(P2_var_on_right_commute_list)
                            sub_terms_list_copy_1[LHS_commutator_P2_index] = f"a2{left_P2_commutator_index}*U2{right_P2_commutator_index}"
                            RHS_commutator_P2_index = sub_terms_list_copy_1.index(P2_var)
                            sub_terms_list_copy_1.pop(RHS_commutator_P2_index)
                            if sub_terms_list_copy_1[0][0] != '-':
                                sub_terms_list_copy_1 = ['-(1/l^2)'] + sub_terms_list_copy_1
                            else:
                                sub_terms_list_copy_1 = ['(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_1))

                            LHS_commutator_P2_index = sub_terms_list_copy_2.index(P2_var_on_right_commute_list)
                            sub_terms_list_copy_2[LHS_commutator_P2_index] = f"a2{right_P2_commutator_index}*U2{left_P2_commutator_index}"
                            RHS_commutator_P2_index = sub_terms_list_copy_2.index(P2_var)
                            sub_terms_list_copy_2.pop(RHS_commutator_P2_index)
                            if sub_terms_list_copy_2[0][0] != '-':
                                sub_terms_list_copy_2 = ['(1/l^2)'] + sub_terms_list_copy_2
                            else:
                                sub_terms_list_copy_2 = ['-(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # lastly swap the orders of P2s in both lists for the next iteration
                            sub_terms_list[P2_var_on_right_sub_list], sub_terms_list[P2_index_sub_list] = (
                                sub_terms_list[P2_index_sub_list],
                                sub_terms_list[P2_var_on_right_sub_list]
                            )
                            P2_to_commute[P2_index_commute_list + 1], P2_to_commute[P2_index_commute_list] = (
                                P2_to_commute[P2_index_commute_list],
                                P2_to_commute[P2_index_commute_list + 1]
                            )
                        # -------------------------------------------------------
                        # Commute second P2
                        # -------------------------------------------------------
                        P22_to_commute = [var for var in sub_terms_list if var.startswith('P2')]
                        while P22_to_commute.index(P22_var) != (len(P22_to_commute) - 1):
                            P22_index_commute_list = P22_to_commute.index(P22_var)
                            P22_index_sub_list = sub_terms_list.index(P22_var)
                            P22_var_on_right_commute_list = P22_to_commute[P22_index_commute_list + 1]
                            P22_var_on_right_sub_list = sub_terms_list.index(P22_var_on_right_commute_list)
                            right_P22_commutator_index = P22_var_on_right_commute_list[-2:]
                            left_P22_commutator_index = P22_var[-2:]
                            P22s_on_right = [
                                var for var in P22_to_commute
                                if P22_to_commute.index(var) > P22_index_commute_list + 1
                            ]
                            if P22_var_on_right_commute_list[-1] == P22_var[-1]:
                                sub_terms_list[P22_var_on_right_sub_list], sub_terms_list[P22_index_sub_list] = (
                                    sub_terms_list[P22_index_sub_list],
                                    sub_terms_list[P22_var_on_right_sub_list]
                                )
                                P22_to_commute[P22_index_commute_list + 1], P22_to_commute[P22_index_commute_list] = (
                                    P22_to_commute[P22_index_commute_list],
                                    P22_to_commute[P22_index_commute_list + 1]
                                )
                                continue
                            for P22_r in P22s_on_right:
                                sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*P2_a*g_cb
                                sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*P2_b*g_ca
                                left_P22_commutator_pos = sub_terms_list_copy_1.index(P22_var_on_right_commute_list)
                                sub_terms_list_copy_1.pop(left_P22_commutator_pos)
                                right_P22_commutator_pos = sub_terms_list_copy_1.index(P22_var)
                                sub_terms_list_copy_1.pop(right_P22_commutator_pos)
                                sub_terms_list_copy_1[sub_terms_list_copy_1.index(P22_r)] = f"P2{left_P22_commutator_index}"
                                if P22_var_on_right_commute_list[-1] == P22_r[-1]:
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['d*(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))
                                else:
                                    for var in sub_terms_list_copy_1:
                                        if var[-1] == P22_r[-1]:
                                            sub_list_var_index = sub_terms_list_copy_1.index(var)
                                            sub_terms_list_copy_1[sub_list_var_index] = f"{var[:-2]}{right_P22_commutator_index}"
                                            break
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))

                                # second term
                                left_P22_commutator_pos = sub_terms_list_copy_2.index(P22_var_on_right_commute_list)
                                sub_terms_list_copy_2.pop(left_P22_commutator_pos)
                                right_P3_commutator_pos = sub_terms_list_copy_2.index(P22_var)
                                sub_terms_list_copy_2.pop(right_P3_commutator_pos)
                                sub_terms_list_copy_2[sub_terms_list_copy_2.index(P22_r)] = f"P2{right_P22_commutator_index}"
                                if P22_var[-1] == P22_r[-1]:
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-d*(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['d*(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                                else:
                                    for var in sub_terms_list_copy_2:
                                        if var[-1] == P22_r[-1]:
                                            sub_list_var_index = sub_terms_list_copy_2.index(var)
                                            sub_terms_list_copy_2[sub_list_var_index] = f"{var[:-2]}{left_P22_commutator_index}"
                                            break
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # now we have to add in the riemann tensor operator that acts on the fields themselves.
                            sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*a1_b*U2_a
                            sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*a1_a*U2_b
                            LHS_commutator_P22_index = sub_terms_list_copy_1.index(P22_var_on_right_commute_list)
                            sub_terms_list_copy_1[LHS_commutator_P22_index] = (
                                f"a2{left_P22_commutator_index}*U2{right_P22_commutator_index}"
                            )
                            RHS_commutator_P22_index = sub_terms_list_copy_1.index(P22_var)
                            sub_terms_list_copy_1.pop(RHS_commutator_P22_index)
                            if sub_terms_list_copy_1[0][0] != '-':
                                sub_terms_list_copy_1 = ['-(1/l^2)'] + sub_terms_list_copy_1
                            else:
                                sub_terms_list_copy_1 = ['(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_1))

                            LHS_commutator_P22_index = sub_terms_list_copy_2.index(P22_var_on_right_commute_list)
                            sub_terms_list_copy_2[LHS_commutator_P22_index] = (
                                f"a2{right_P22_commutator_index}*U2{left_P22_commutator_index}"
                            )
                            RHS_commutator_P22_index = sub_terms_list_copy_2.index(P22_var)
                            sub_terms_list_copy_2.pop(RHS_commutator_P22_index)
                            if sub_terms_list_copy_2[0][0] != '-':
                                sub_terms_list_copy_2 = ['(1/l^2)'] + sub_terms_list_copy_2
                            else:
                                sub_terms_list_copy_2 = ['-(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # lastly swap the orders of P2s in both lists for the next iteration
                            sub_terms_list[P22_var_on_right_sub_list], sub_terms_list[P22_index_sub_list] = (
                                sub_terms_list[P22_index_sub_list],
                                sub_terms_list[P22_var_on_right_sub_list]
                            )
                            P22_to_commute[P22_index_commute_list + 1], P22_to_commute[P22_index_commute_list] = (
                                P22_to_commute[P22_index_commute_list],
                                P22_to_commute[P22_index_commute_list + 1]
                            )
                        # only want to move one instance - break out of two loops and go to next term
                        break

        if B2_present:
            # now that sub_terms_list has the variables in the desired order,
            # we can replace B2 with m2 
            sub_terms_list_B2 = sub_terms_list.copy()
            # can just use P2_var and P22_var - in first list, delete the P2_var and replace the P22_var with the box operator (same index)
            P2_var_index = sub_terms_list.index(P2_var_stored)
            # delete P2_var in list
            sub_terms_list_B2.pop(P2_var_index)
            # find P22_var index
            P22_var_index_B2 = sub_terms_list_B2.index(P22_var)
            # and also delete
            sub_terms_list_B2.pop(P22_var_index_B2)
            # have to include extra on shell constant m2
            if sub_terms_list_B2[0][0] == '-':
                sub_terms_list_B2 = ['-m2'] + [sub_terms_list_B2[0][1:]] + sub_terms_list_B2[1:]
            else:
                sub_terms_list_B2 = ['m2'] + sub_terms_list_B2
            
            # then we can construct the terms by joining with '*'
            expanded_terms.append('*'.join(sub_terms_list_B2))
        else:
            # just reconstruct the original term
            expanded_terms.append('*'.join(sub_terms_list))

    return pull_minus_signs_to_front(' + '.join(expanded_terms))

def pull_B3s_right(equation):
    # pulls the leftmost B3 to the right to impose on shell conditions
    # break the equation down into individual terms by top level '+' symbols
    equation = replace_minuses(equation)
    terms = top_level_split(equation)
    expanded_terms = []
    for term in terms:
        # break into individual variables by * signs
        # this will be the main list which will be rearranged as we go
        sub_terms_list = term.split('*')
        # check over sub_terms_list for P3 variables
        B3_present = False
        for P3_var in reversed(sub_terms_list):
            if B3_present == True:
                break
            elif P3_var.startswith('P3'):
                P3_var_stored = P3_var
                # store the letter index and position index of the P3
                P3_lett_index = P3_var[-1]
                P3_pos_index = sub_terms_list.index(P3_var)
                for P33_var in reversed(sub_terms_list):
                    if P33_var.startswith('P3') and P33_var[-1] == P3_lett_index and P33_var != P3_var:
                        B3_present = True
                        # store the P33 index (letter index should be same as P3 if relevant)
                        P33_pos_index = sub_terms_list.index(P33_var)

                        # -------------------------------------------------------
                        # Now we want to move all U3, Z2, Z1 left of the leftmost P
                        # before we can commute P3s
                        # -------------------------------------------------------
                        U3_Z2_Z1_list = [var for var in sub_terms_list if var.startswith('U1') or var.startswith('U2') or var.startswith('U3') or var.startswith('Z2') or var.startswith('Z1') or var.startswith('Z3')]
                        # Create a list with those items removed
                        remaining_list = [
                            var for var in sub_terms_list
                            if var not in U3_Z2_Z1_list
                        ]

                        # Find the position of the leftmost P in the remaining_list
                        for var in remaining_list:
                            if var.startswith('P'):
                                leftmost_P = var
                                leftmost_P_pos = remaining_list.index(var)
                                break

                        # Insert all U3/Z2/Z1 vars to the left of that leftmost P
                        for U3_Z2_Z1 in U3_Z2_Z1_list:
                            remaining_list.insert(leftmost_P_pos, U3_Z2_Z1)

                        # Reassign sub_terms_list
                        sub_terms_list = remaining_list

                        # -------------------------------------------------------
                        # Commute P3s
                        # -------------------------------------------------------
                        P3_to_commute = [var for var in sub_terms_list if var.startswith('P3')]
                        while P3_to_commute.index(P3_var) != (len(P3_to_commute) - 1):
                            # find P3 index in the to_commute list
                            P3_index_commute_list = P3_to_commute.index(P3_var)
                            # find P3 index in the sub_terms_list
                            P3_index_sub_list = sub_terms_list.index(P3_var)
                            # store the P3 variable on the next right to be commuted through
                            P3_var_on_right_commute_list = P3_to_commute[P3_index_commute_list + 1]
                            # find this variable's position in the sub_terms_list
                            P3_var_on_right_sub_list = sub_terms_list.index(P3_var_on_right_commute_list)
                            right_P3_commutator_index = P3_var_on_right_commute_list[-2:]
                            left_P3_commutator_index = P3_var[-2:]
                            # create terms arising from the riemann tensor contracting with P3s on the right
                            P3s_on_right = [
                                var for var in P3_to_commute
                                if P3_to_commute.index(var) > P3_index_commute_list + 1
                            ]
                            # if letter index on P3 being commuted is the same, just swap
                            if P3_var_on_right_commute_list[-1] == P3_var[-1]:
                                sub_terms_list[P3_var_on_right_sub_list], sub_terms_list[P3_index_sub_list] = (
                                    sub_terms_list[P3_index_sub_list],
                                    sub_terms_list[P3_var_on_right_sub_list]
                                )
                                P3_to_commute[P3_index_commute_list + 1], P3_to_commute[P3_index_commute_list] = (
                                    P3_to_commute[P3_index_commute_list],
                                    P3_to_commute[P3_index_commute_list + 1]
                                )
                                continue
                            for P3_r in P3s_on_right:
                                # create copy to manipulate terms
                                sub_terms_list_copy_1 = sub_terms_list.copy()   # for making 1/l^2*P3_a*g_cb
                                sub_terms_list_copy_2 = sub_terms_list.copy()   # for making -1/l^2*P3_b*g_ca
                                # First term
                                sub_terms_list_copy_1[sub_terms_list_copy_1.index(P3_r)] = f"P3{left_P3_commutator_index}"
                                left_P3_commutator_pos = sub_terms_list_copy_1.index(P3_var_on_right_commute_list)
                                sub_terms_list_copy_1.pop(left_P3_commutator_pos)
                                right_P3_commutator_pos = sub_terms_list_copy_1.index(P3_var)
                                sub_terms_list_copy_1.pop(right_P3_commutator_pos)

                                if P3_var_on_right_commute_list[-1] == P3_r[-1]:
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['d*(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))
                                else:
                                    for var in sub_terms_list_copy_1:
                                        if var[-1] == P3_r[-1]:
                                            sub_list_var_index = sub_terms_list_copy_1.index(var)
                                            sub_terms_list_copy_1[sub_list_var_index] = f"{var[:-2]}{right_P3_commutator_index}"
                                            break
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))

                                # Second term
                                left_P3_commutator_pos = sub_terms_list_copy_2.index(P3_var_on_right_commute_list)
                                sub_terms_list_copy_2.pop(left_P3_commutator_pos)
                                right_P3_commutator_pos = sub_terms_list_copy_2.index(P3_var)
                                sub_terms_list_copy_2.pop(right_P3_commutator_pos)
                                sub_terms_list_copy_2[sub_terms_list_copy_2.index(P3_r)] = f"P3{right_P3_commutator_index}"

                                if P3_var[-1] == P3_r[-1]:
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['-d*(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))
                                else:
                                    for var in sub_terms_list_copy_2:
                                        if var[-1] == P3_r[-1]:
                                            sub_list_var_index = sub_terms_list_copy_2.index(var)
                                            sub_terms_list_copy_2[sub_list_var_index] = f"{var[:-2]}{left_P3_commutator_index}"
                                            break
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # Riemann tensor operator on fields
                            sub_terms_list_copy_1 = sub_terms_list.copy()
                            sub_terms_list_copy_2 = sub_terms_list.copy()
                            LHS_commutator_P3_index = sub_terms_list_copy_1.index(P3_var_on_right_commute_list)
                            sub_terms_list_copy_1[LHS_commutator_P3_index] = f"a3{left_P3_commutator_index}*U3{right_P3_commutator_index}"
                            RHS_commutator_P3_index = sub_terms_list_copy_1.index(P3_var)
                            sub_terms_list_copy_1.pop(RHS_commutator_P3_index)
                            if sub_terms_list_copy_1[0][0] != '-':
                                sub_terms_list_copy_1 = ['-(1/l^2)'] + sub_terms_list_copy_1
                            else:
                                sub_terms_list_copy_1 = ['(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_1))

                            LHS_commutator_P3_index = sub_terms_list_copy_2.index(P3_var_on_right_commute_list)
                            sub_terms_list_copy_2[LHS_commutator_P3_index] = f"a3{right_P3_commutator_index}*U3{left_P3_commutator_index}"
                            RHS_commutator_P3_index = sub_terms_list_copy_2.index(P3_var)
                            sub_terms_list_copy_2.pop(RHS_commutator_P3_index)
                            if sub_terms_list_copy_2[0][0] != '-':
                                sub_terms_list_copy_2 = ['(1/l^2)'] + sub_terms_list_copy_2
                            else:
                                sub_terms_list_copy_2 = ['-(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # Swap P3s
                            sub_terms_list[P3_var_on_right_sub_list], sub_terms_list[P3_index_sub_list] = (
                                sub_terms_list[P3_index_sub_list],
                                sub_terms_list[P3_var_on_right_sub_list]
                            )
                            P3_to_commute[P3_index_commute_list + 1], P3_to_commute[P3_index_commute_list] = (
                                P3_to_commute[P3_index_commute_list],
                                P3_to_commute[P3_index_commute_list + 1]
                            )
                        # -------------------------------------------------------
                        # Commute second P3
                        # -------------------------------------------------------
                        P33_to_commute = [var for var in sub_terms_list if var.startswith('P3')]
                        while P33_to_commute.index(P33_var) != (len(P33_to_commute) - 1):
                            P33_index_commute_list = P33_to_commute.index(P33_var)
                            P33_index_sub_list = sub_terms_list.index(P33_var)
                            P33_var_on_right_commute_list = P33_to_commute[P33_index_commute_list + 1]
                            P33_var_on_right_sub_list = sub_terms_list.index(P33_var_on_right_commute_list)
                            right_P33_commutator_index = P33_var_on_right_commute_list[-2:]
                            left_P33_commutator_index = P33_var[-2:]
                            P33s_on_right = [
                                var for var in P33_to_commute
                                if P33_to_commute.index(var) > P33_index_commute_list + 1
                            ]
                            if P33_var_on_right_commute_list[-1] == P33_var[-1]:
                                sub_terms_list[P33_var_on_right_sub_list], sub_terms_list[P33_index_sub_list] = (
                                    sub_terms_list[P33_index_sub_list],
                                    sub_terms_list[P33_var_on_right_sub_list]
                                )
                                P33_to_commute[P33_index_commute_list + 1], P33_to_commute[P33_index_commute_list] = (
                                    P33_to_commute[P33_index_commute_list],
                                    P33_to_commute[P33_index_commute_list + 1]
                                )
                                continue
                            for P33_r in P33s_on_right:
                                sub_terms_list_copy_1 = sub_terms_list.copy()
                                sub_terms_list_copy_2 = sub_terms_list.copy()
                                left_P33_commutator_pos = sub_terms_list_copy_1.index(P33_var_on_right_commute_list)
                                sub_terms_list_copy_1.pop(left_P33_commutator_pos)
                                right_P33_commutator_pos = sub_terms_list_copy_1.index(P33_var)
                                sub_terms_list_copy_1.pop(right_P33_commutator_pos)
                                sub_terms_list_copy_1[sub_terms_list_copy_1.index(P33_r)] = f"P3{left_P33_commutator_index}"
                                if P33_var_on_right_commute_list[-1] == P33_r[-1]:
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['d*(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-d*(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))
                                else:
                                    for var in sub_terms_list_copy_1:
                                        if var[-1] == P33_r[-1]:
                                            i_var = sub_terms_list_copy_1.index(var)
                                            sub_terms_list_copy_1[i_var] = f"{var[:-2]}{right_P33_commutator_index}"
                                            break
                                    if sub_terms_list_copy_1[0][0] != '-':
                                        sub_terms_list_copy_1 = ['(1/l^2)'] + sub_terms_list_copy_1
                                    else:
                                        sub_terms_list_copy_1 = ['-(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_1))

                                # second term
                                left_P33_commutator_pos = sub_terms_list_copy_2.index(P33_var_on_right_commute_list)
                                sub_terms_list_copy_2.pop(left_P33_commutator_pos)
                                right_P3_commutator_pos = sub_terms_list_copy_2.index(P33_var)
                                sub_terms_list_copy_2.pop(right_P3_commutator_pos)
                                sub_terms_list_copy_2[sub_terms_list_copy_2.index(P33_r)] = f"P3{right_P33_commutator_index}"
                                if P33_var[-1] == P33_r[-1]:
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-d*(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['d*(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))
                                else:
                                    for var in sub_terms_list_copy_2:
                                        if var[-1] == P33_r[-1]:
                                            i_var = sub_terms_list_copy_2.index(var)
                                            sub_terms_list_copy_2[i_var] = f"{var[:-2]}{left_P33_commutator_index}"
                                            break
                                    if sub_terms_list_copy_2[0][0] != '-':
                                        sub_terms_list_copy_2 = ['-(1/l^2)'] + sub_terms_list_copy_2
                                    else:
                                        sub_terms_list_copy_2 = ['(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                                    expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # now we have to add in the riemann tensor operator that acts on the fields themselves.
                            sub_terms_list_copy_1 = sub_terms_list.copy()
                            sub_terms_list_copy_2 = sub_terms_list.copy()
                            LHS_commutator_P33_index = sub_terms_list_copy_1.index(P33_var_on_right_commute_list)
                            sub_terms_list_copy_1[LHS_commutator_P33_index] = (
                                f"a3{left_P33_commutator_index}*U3{right_P33_commutator_index}"
                            )
                            RHS_commutator_P33_index = sub_terms_list_copy_1.index(P33_var)
                            sub_terms_list_copy_1.pop(RHS_commutator_P33_index)
                            if sub_terms_list_copy_1[0][0] != '-':
                                sub_terms_list_copy_1 = ['-(1/l^2)'] + sub_terms_list_copy_1
                            else:
                                sub_terms_list_copy_1 = ['(1/l^2)'] + [sub_terms_list_copy_1[0][1:]] + sub_terms_list_copy_1[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_1))

                            LHS_commutator_P33_index = sub_terms_list_copy_2.index(P33_var_on_right_commute_list)
                            sub_terms_list_copy_2[LHS_commutator_P33_index] = (
                                f"a3{right_P33_commutator_index}*U3{left_P33_commutator_index}"
                            )
                            RHS_commutator_P33_index = sub_terms_list_copy_2.index(P33_var)
                            sub_terms_list_copy_2.pop(RHS_commutator_P33_index)
                            if sub_terms_list_copy_2[0][0] != '-':
                                sub_terms_list_copy_2 = ['(1/l^2)'] + sub_terms_list_copy_2
                            else:
                                sub_terms_list_copy_2 = ['-(1/l^2)'] + [sub_terms_list_copy_2[0][1:]] + sub_terms_list_copy_2[1:]
                            expanded_terms.append('*'.join(sub_terms_list_copy_2))

                            # Swap the orders of P3s
                            sub_terms_list[P33_var_on_right_sub_list], sub_terms_list[P33_index_sub_list] = (
                                sub_terms_list[P33_index_sub_list],
                                sub_terms_list[P33_var_on_right_sub_list]
                            )
                            P33_to_commute[P33_index_commute_list + 1], P33_to_commute[P33_index_commute_list] = (
                                P33_to_commute[P33_index_commute_list],
                                P33_to_commute[P33_index_commute_list + 1]
                            )
                        # only want to move one instance - break out of two loops
                        break

        if B3_present:
            # now that sub_terms_list has the variables in the desired order,
            # we can replace B3 with m3 
            sub_terms_list_B3 = sub_terms_list.copy()
            # can just use P3_var and P33_var - in first list, delete the P3_var and replace the P33_var with the box operator (same index)
            P3_var_index = sub_terms_list.index(P3_var_stored)
            # delete P3_var in list
            sub_terms_list_B3.pop(P3_var_index)
            # find P33_var index
            P33_var_index_B3 = sub_terms_list_B3.index(P33_var)
            # and also delete
            sub_terms_list_B3.pop(P33_var_index_B3)
            # have to include extra on shell constant m3
            if sub_terms_list_B3[0][0] == '-':
                sub_terms_list_B3 = ['-m3'] + [sub_terms_list_B3[0][1:]] + sub_terms_list_B3[1:]
            else:
                sub_terms_list_B3 = ['m3'] + sub_terms_list_B3
            
            # then we can construct the terms by joining with '*'
            expanded_terms.append('*'.join(sub_terms_list_B3))
        else:
            # just reconstruct the original term
            expanded_terms.append('*'.join(sub_terms_list))

    return pull_minus_signs_to_front(' + '.join(expanded_terms))

def pull_all_Bs_right(equation):
    return move_all_a_left(pull_B1s_right(move_all_a_left(pull_B2s_right(move_all_a_left(pull_B3s_right(equation))))))

def impose_y_rule(equation):
    terms = equation.split(" + ")
    processed_terms = []
    for term in terms:
        # Impose rules for Y1 operator
        term = re.sub(r"U1\_([a-z])\*P2\^\1", "Y1", term)
        term = re.sub(r"P2\_([a-z])\*U1\^\1", "Y1", term)
        term = re.sub(r"U1\^([a-z])\*P2\_\1", "Y1", term)
        term = re.sub(r"P2\^([a-z])\*U1\_\1", "Y1", term)

        # Impose rules for Y2 operator
        term = re.sub(r"U2\_([a-z])\*P3\^\1", "Y2", term)
        term = re.sub(r"P3\_([a-z])\*U2\^\1", "Y2", term)
        term = re.sub(r"U2\^([a-z])\*P3\_\1", "Y2", term)
        term = re.sub(r"P3\^([a-z])\*U2\_\1", "Y2", term)

        # Impose rules for Y3 operators
        term = re.sub(r"U3\_([a-z])\*P1\^\1", "Y3", term)
        term = re.sub(r"P1\_([a-z])\*U3\^\1", "Y3", term)
        term = re.sub(r"U3\^([a-z])\*P1\_\1", "Y3", term)
        term = re.sub(r"P1\^([a-z])\*U3\_\1", "Y3", term)

        
        processed_terms.append(term)

    return " + ".join(processed_terms)

def impose_z_rule(equation):
    """
    Processes an equation string by replacing specific contractions of U1, U2, and U3 
    with Z3, Z2, and Z1 respectively, while iteratively applying the rules.
    
    Parameters:
        equation (str): The input equation as a single string.

    Returns:
        str: The processed equation.
    """
    # Split terms by '+' (keeping terms distinct)
    terms_list = equation.split(' + ')

    def replace_contractions(term):
        prev_term = None
        while term != prev_term:  # Keep processing until no further changes
            prev_term = term
            # Replace contractions of U1 and U2 with Z3
            term = re.sub(r"U1_([a-z])\*U2\^\1", "Z3", term)
            term = re.sub(r"U2_([a-z])\*U1\^\1", "Z3", term)
            term = re.sub(r"U1\^([a-z])\*U2\_\1", "Z3", term)
            term = re.sub(r"U2\^([a-z])\*U1\_\1", "Z3", term)

            # Replace contractions of U1 and U3 with Z2
            term = re.sub(r"U1_([a-z])\*U3\^\1", "Z2", term)
            term = re.sub(r"U3_([a-z])\*U1\^\1", "Z2", term)
            term = re.sub(r"U1\^([a-z])\*U3\_\1", "Z2", term)
            term = re.sub(r"U3\^([a-z])\*U1\_\1", "Z2", term)

            # Replace contractions of U2 and U3 with Z1
            term = re.sub(r"U2_([a-z])\*U3\^\1", "Z1", term)
            term = re.sub(r"U3_([a-z])\*U2\^\1", "Z1", term)
            term = re.sub(r"U2\^([a-z])\*U3\_\1", "Z1", term)
            term = re.sub(r"U3\^([a-z])\*U2\_\1", "Z1", term)
        return term

    # Process all terms
    processed_terms = [replace_contractions(term.strip()) for term in terms_list]

    # Rebuild the equation
    return " + ".join(processed_terms)

def impose_dimension(equation, d=3):
    # Substitute in a particular dimension (default d=3) into equations,
    # but do not replace 'd' if it is immediately preceded by '_' or '^'.
    # We use a negative lookbehind: (?<![_\^])d matches a "d" that is not preceded by "_" or "^".
    terms = top_level_split(equation)
    new_terms = []
    pattern = r'(?<![_\^])d'
    for term in terms:
        new_term = re.sub(pattern, str(d), term)
        new_terms.append(new_term)
    return ' + '.join(new_terms)

def combine_same_terms(equation):
    """
    Combines like terms by parsing numeric factors as Fractions,
    summing them if leftover factors match exactly, and returning
    a final expression with fractional coefficients displayed as (3/2) etc.

    Key points:
      - We do NOT expand (s3-1) or unify exponents.  We treat them as leftover text.
      - If net coefficient is a fraction like 1/2, we show (1/2).
      - If net coefficient is integer, we show e.g. 2 or -3 or 1.
    """

    # 1) We assume top_level_split(equation) => list of raw terms (+/-).
    terms = top_level_split(equation)

    # Patterns to detect purely numeric factors:
    #   integer/decimal or strictly numeric fraction in parentheses
    purely_digit_pattern = re.compile(r'^[+\-]?\d+(?:\.\d+)?$')
    fraction_pattern = re.compile(r'''
        ^
        \(\s*
        [+\-]?\d+(?:\.\d+)?    # optional sign, digits, optional decimal
        \s*/\s*
        [+\-]?\d+(?:\.\d+)?    # optional sign, digits, optional decimal
        \s*\)
        $
    ''', re.VERBOSE)

    def is_numeric_factor(fac):
        """Return True if 'fac' is recognized as numeric: integer/float or (num/num) with no letters."""
        test = fac.strip()
        if purely_digit_pattern.match(test):
            return True
        if fraction_pattern.match(test):
            return True
        return False

    def parse_as_fraction(num_str):
        """
        Convert '3','-2.5' or '(3/4)' to Fraction.
        If (x/y), remove parentheses and parse x/y. If invalid => 0.
        """
        s = num_str.strip()
        if s.startswith('(') and s.endswith(')') and '/' in s:
            inner = s[1:-1].strip()  # e.g. "3/4"
            return Fraction(inner)
        else:
            try:
                return Fraction(s)
            except ValueError:
                return Fraction(0)

    # For leftover factors, we reorder so non-Z/Y first => alpha, then Z or Y => alpha
    def factor_key(f):
        if f.startswith('Z') or f.startswith('Y'):
            return (1, f)
        else:
            return (0, f)

    def canonicalize_term(raw_term):
        """
        1) parse leading sign => ±1
        2) split by '*'
        3) multiply numeric => net fraction
        4) leftover => reorder => leftover_str
        5) return (Fraction_coefficient, leftover_str)
        """
        s = raw_term.strip()
        sign = 1
        if s.startswith('-'):
            sign = -1
            s = s[1:].strip()
        elif s.startswith('+'):
            s = s[1:].strip()

        pieces = [x.strip() for x in s.split('*') if x.strip()]
        coeff = Fraction(sign, 1)
        leftover = []

        for p in pieces:
            if is_numeric_factor(p):
                coeff *= parse_as_fraction(p)
            else:
                leftover.append(p)

        leftover.sort(key=factor_key)
        leftover_str = '*'.join(leftover)
        return (coeff, leftover_str)

    # (A) Collect net fraction for each leftover_str
    store = defaultdict(Fraction)
    for t in terms:
        c, l_str = canonicalize_term(t)
        store[l_str] += c

    # remove zeros
    for k in list(store.keys()):
        if store[k] == 0:
            del store[k]

    if not store:
        return '0'

    # sort by leftover_str for stable output
    sorted_items = sorted(store.items(), key=lambda x: x[0])

    # (B) Build final expression
    final_terms = []
    for leftover_str, val in sorted_items:
        if val == 0:
            continue
        sign_char = ''
        abs_val = val
        if val < 0:
            sign_char = '-'
            abs_val = -val

        # if leftover_str is empty => purely numeric product
        # if abs_val=1 => omit numeric if leftover
        if abs_val == 1:
            if leftover_str:
                piece = leftover_str
            else:
                piece = '1'
        else:
            # Convert fraction to string
            abs_str = str(abs_val)  # e.g. "1/2" or "3/4" or "2" or "5"
            # If denominator != 1 => it's a fraction => wrap in parentheses
            if abs_val.denominator != 1:
                abs_str = f"({abs_str})"
            if leftover_str:
                piece = abs_str + '*' + leftover_str
            else:
                piece = abs_str

        if sign_char == '-':
            final_terms.append('-' + piece)
        else:
            final_terms.append(piece)

    if not final_terms:
        return '0'

    # (C) Join with ' + ', turning e.g. "stuff + -stuff2" => "stuff - stuff2"
    out_list = []
    for i, term in enumerate(final_terms):
        t = term.strip()
        if i == 0:
            out_list.append(t)  # no leading plus if positive
        else:
            out_list.append('+ ' + t)

    return ' '.join(out_list)

def expand_a_bracket(equation):
    # try to implement more general method for expanding brackets
    terms = top_level_split(equation)
    new_terms = []
    for term in terms:
        term_is_neg = term[0] == '-'
        if term_is_neg:     
            sub_terms_list = top_level_split(term[1:], delimiter='*')
        else:
            sub_terms_list = top_level_split(term, delimiter='*')
        # want to see if there are terms containing brackets with more sub-terms
        for var in sub_terms_list:
            if var[0] == '(' and var[-1] == ')':
                split_var = top_level_split_sign(var[1:-1])
            else:
                split_var = top_level_split_sign(var)
            split_var_no_pluses = []
            for subvar in split_var:
                if subvar[0] == '+':
                    split_var_no_pluses.append(subvar[1:])
                else:
                    split_var_no_pluses.append(subvar)
            no_sub_brackets = True
            if len(split_var_no_pluses) > 1:
                # if we find a bracketed term with sub terms, we want to generate new terms, one for each of all the other terms multiplied by each
                # variable in the bracketed sub term (expanding brackets) 
                for sub_var in split_var_no_pluses:
                    sub_terms_list_copy = sub_terms_list.copy()
                    # in the copy, replace the bracket term with one term from the bracket (repeated for each subterm)
                    orig_bracket_pos = sub_terms_list_copy.index(var)
                    sub_terms_list_copy[orig_bracket_pos] = sub_var
                    if term_is_neg:
                        new_term = '-1*' + '*'.join(sub_terms_list_copy)
                    else:
                        new_term = '*'.join(sub_terms_list_copy)
                    new_terms.append(new_term)
                no_sub_brackets = False
                # only do it for one bracket at a time
                break
        if no_sub_brackets == True:
            if term_is_neg:
                new_term = '-1*' + '*'.join(sub_terms_list)
            else:
                new_term = '*'.join(sub_terms_list)
            new_terms.append(new_term)
    return pull_minus_signs_to_front(combine_numerical_factors(pull_minus_signs_to_front(' + '.join(new_terms))))

def combine_numerical_factors(equation):
    """
    For each top‐level term (as returned by top_level_split), split it on '*' and identify numeric factors.
    A numeric factor is defined as either:
      1) A standalone integer/decimal (e.g. "2", "-3.5")
      2) A parenthesized fraction that is entirely numeric (e.g. "(3/4)"); now also allowing an optional sign before the parentheses,
         e.g. "-(1/2)" or "+(3/2)".
    Factors that include any letters (e.g. "s2" or "(s3-1)") are not treated as numeric.
    
    For each term, the function prints the term together with the list of numerical factors it found.
    """
    # Assume top_level_split is defined elsewhere.
    terms = top_level_split(equation)
    
    # Pattern for integer or decimal numbers (standalone)
    # Matches: 2, -3.5, +4, etc.
    integer_decimal_pattern = re.compile(
        r'(?<![\w])'            # not preceded by a letter, digit, or underscore
        r'[+\-]?\d+(?:\.\d+)?'   # digits with optional sign and optional decimal
        r'(?![\w])'             # not followed by a letter, digit, or underscore
    )
    
    # Pattern for a parenthesized fraction. Now we allow an optional sign BEFORE the parentheses.
    # It matches strings like: "(3/4)", "-(1/2)", "+(3/2)", etc.
    fraction_pattern = re.compile(
        r'(?<![\w])'                     # not preceded by letter/digit/underscore
        r'[+\-]?'                        # optional sign
        r'\(\s*[+\-]?\d+(?:\.\d+)?\s*/\s*[+\-]?\d+(?:\.\d+)?\s*\)'  # a fraction in parentheses (only digits inside)
        r'(?![\w])'                      # not followed by letter/digit/underscore
    )
    
    def is_numeric_factor(fac):
        fac = fac.strip()
        # use fullmatch so that the entire string must match our pattern
        if integer_decimal_pattern.fullmatch(fac):
            return True
        if fraction_pattern.fullmatch(fac):
            return True
        return False
    
    new_terms = []
    for term in terms:
        # Split the term by '*' into factors.
        factors = [f.strip() for f in top_level_split(term, delimiter='*') if f.strip()]
        # Pick out those factors that are purely numeric
        numeric_factors = [fac for fac in factors if is_numeric_factor(fac)]
        non_numeric_factors = [f for f in factors if f not in numeric_factors]
        result = Fraction(1, 1)
        for f in numeric_factors:
            s = f.strip()
            # Use regex to remove outer parentheses, if present.
            # This pattern looks for an optional leading sign, then a '(',
            # some content (captured), and then a ')', and replaces it with sign+content.
            s = re.sub(r'^([+\-]?)\(\s*(.*?)\s*\)$', r'\1\2', s)
            result *= Fraction(s)
        result = str(result)
        if '/' in result:
            if result[0] == '-':
                result = '-(' + result[1:] + ')'
            else:
                result = '(' + result + ')'
        term = '*'.join(non_numeric_factors)
        if result == '1':
            if term == '':
                new_terms.append('1')
            else:
                new_terms.append(term)
        elif result == '-1':
            if term == '':
                new_terms.append('-1')
            else:
                new_terms.append('-' + term)
        else:
            new_terms.append(result + '*' + term)
    return ' + '.join(new_terms) 

def fully_expand_equation(equation):
    while equation != expand_a_bracket(equation):
        equation = expand_a_bracket(equation)
    return combine_numerical_factors(equation)

def substitute_n1_n2_n3(equation_str):
    """
    Takes an expanded equation string, e.g.:
      '(1/l^2)*n3*s3*Z1^(n1+2) + -(1/l^2)*n3*Z2^(n2+1)*Z3^(n3+1) + ...'
    and substitutes:
      n1 -> (1/2)*(s2+s3-s1-2)
      n2 -> (1/2)*(s3+s1-s2-2)
      n3 -> (1/2)*(s1+s2-s3-2)
    ONLY in factors that are NOT exponents of the form 'Z...^(...)' or 'Y...^(...)'.
    Additionally, strips any directly surrounding brackets around n1, n2, n3 before substitution.
    
    Summands are assumed to be separated by ' + ', and negative terms
    might appear as '-stuff' or '+ -stuff'. The function:
      1) Splits by ' + '.
      2) Checks for a leading '-'.
      3) Splits each chunk by '*'.
      4) Strips surrounding brackets from n1, n2, n3 if present.
      5) If factor doesn't match r'^[ZY]\\d+\\^', performs substitution.
      6) Rejoins factors with '*'.
      7) Reattaches '-' if originally negative.
      8) Rejoins all summands with ' + ' again.
    """

    # The strings to substitute
    sub_n1 = '(1/2)*(s2+s3-s1-2)'
    sub_n2 = '(1/2)*(s3+s1-s2-2)'
    sub_n3 = '(1/2)*(s1+s2-s3-2)'

    # 1) Split the entire expression by " + "
    #    We assume your final expansions are separated by " + " (some might start with '-').
    #    e.g.: "stuff + -otherStuff + moreStuff"
    summands = equation_str.strip().split(' + ')

    new_summands = []

    for summand in summands:
        s = summand.strip()
        is_negative = False
        if s.startswith('-'):
            # Then the real factor is s[1:]
            is_negative = True
            s = s[1:].strip()

        # 2) Split by '*'
        factors = [f.strip() for f in s.split('*') if f.strip()]

        new_factors = []
        for fac in factors:
            # 3) Strip surrounding brackets if the factor is exactly (n1), (n2), or (n3)
            fac = re.sub(r'^\((n1|n2|n3)\)$', r'\1', fac)

            # 4) Skip substitution if factor looks like Z..^(...) or Y..^(...)
            if re.match(r'^[ZY]\d+\^', fac):
                new_factors.append(fac)
            else:
                # 5) Perform substitution for n1, n2, n3 with word boundaries
                tmp = fac
                tmp = re.sub(r'\bn1\b', sub_n1, tmp)
                tmp = re.sub(r'\bn2\b', sub_n2, tmp)
                tmp = re.sub(r'\bn3\b', sub_n3, tmp)
                new_factors.append(tmp)

        # 6) Rejoin factors with '*'
        joined_factors = '*'.join(new_factors)
        if is_negative:
            joined_factors = '-' + joined_factors
        new_summands.append(joined_factors)

    # 7) Rejoin all summands with " + "
    return ' + '.join(new_summands)

def substitute_p1_p2_p3(equation_str):
    """
    Takes an expanded equation string, e.g.:
      '(1/l^2)*p3*s3*Z1^(p1+2) + -(1/l^2)*p3*Z2^(p2+1)*Z3^(p3+1) + ...'
    and substitutes:
      p1 -> (1/2)*(s2+s3-s1-1)
      p2 -> (1/2)*(s3+s1-s2-1)
      p3 -> (1/2)*(s1+s2-s3-1)
    ONLY in factors that are NOT exponents of the form 'Z...^(...)' or 'Y...^(...)'.
    Additionally, strips any directly surrounding brackets around n1, n2, n3 before substitution.
    
    Summands are assumed to be separated by ' + ', and negative terms
    might appear as '-stuff' or '+ -stuff'. The function:
      1) Splits by ' + '.
      2) Checks for a leading '-'.
      3) Splits each chunk by '*'.
      4) Strips surrounding brackets from p1, p2, p3 if present.
      5) If factor doesn't match r'^[ZY]\\d+\\^', performs substitution.
      6) Rejoins factors with '*'.
      7) Reattaches '-' if originally negative.
      8) Rejoins all summands with ' + ' again.
    """

    # The strings to substitute
    sub_p1 = '(1/2)*(s2+s3-s1-1)'
    sub_p2 = '(1/2)*(s3+s1-s2-1)'
    sub_p3 = '(1/2)*(s1+s2-s3-1)'

    # 1) Split the entire expression by " + "
    #    We assume your final expansions are separated by " + " (some might start with '-').
    #    e.g.: "stuff + -otherStuff + moreStuff"
    summands = equation_str.strip().split(' + ')

    new_summands = []

    for summand in summands:
        s = summand.strip()
        is_negative = False
        if s.startswith('-'):
            # Then the real factor is s[1:]
            is_negative = True
            s = s[1:].strip()

        # 2) Split by '*'
        factors = [f.strip() for f in s.split('*') if f.strip()]

        new_factors = []
        for fac in factors:
            # 3) Strip surrounding brackets if the factor is exactly (p1), (p2), or (p3)
            fac = re.sub(r'^\((p1|p2|p3)\)$', r'\1', fac)

            # 4) Skip substitution if factor looks like Z..^(...) or Y..^(...)
            if re.match(r'^[ZY]\d+\^', fac):
                new_factors.append(fac)
            else:
                # 5) Perform substitution for p1, p2, p3 with word boundaries
                tmp = fac
                tmp = re.sub(r'\bp1\b', sub_p1, tmp)
                tmp = re.sub(r'\bp2\b', sub_p2, tmp)
                tmp = re.sub(r'\bp3\b', sub_p3, tmp)
                new_factors.append(tmp)

        # 6) Rejoin factors with '*'
        joined_factors = '*'.join(new_factors)
        if is_negative:
            joined_factors = '-' + joined_factors
        new_summands.append(joined_factors)

    # 7) Rejoin all summands with " + "
    return ' + '.join(new_summands)

def substitute_m_variables(equation, field_is_gauge_parameter=[1]):
    """
    Given an equation string, substitutes:
       m1 -> (1/l^2)*s1*(s1-1) for gauge parameters,
       m1 -> (1/l^2)*s1*(s1-3) for non‑gauge parameters,
       similarly for m2 and m3.
    Only standalone occurrences (using word boundaries) are replaced.
    """
    import re
    eq = equation  # start with the original

    # For each field number that is a gauge parameter, substitute with (1/l^2)*s{field}*(s{field}-1)
    for field_number in field_is_gauge_parameter:
        eq = re.sub(rf'\bm{field_number}\b', 
                    rf'(1/l^2)*s{field_number}*(s{field_number}-1)', 
                    eq)

    # For each field number not marked as gauge parameter, substitute with (1/l^2)*s{field}*(s{field}-3)
    for field_number in [1, 2, 3]:
        if field_number not in field_is_gauge_parameter:
            eq = re.sub(rf'\bm{field_number}\b', 
                        rf'(1/l^2)*s{field_number}*(s{field_number}-3)', 
                        eq)
    return eq

def complete_sub_and_expansion(equation):
    equation = fully_expand_equation(fully_expand_equation(fully_expand_equation(equation)))
    equation = combine_same_terms(equation)
    equation = substitute_m_variables(equation)
    equation = fully_expand_equation(equation)
    equation = combine_same_terms(equation)
    equation = substitute_n1_n2_n3(equation)
    equation = fully_expand_equation(equation)
    equation = substitute_p1_p2_p3(equation)
    equation = fully_expand_equation(equation)
    equation = combine_same_terms(equation)
    equation = impose_dimension(equation)
    equation = fully_expand_equation(equation)
    equation = combine_same_terms(combine_powers_in_terms(equation))
    return equation

def antisymmetrise_down_indices(equation):
    """
    Given an equation string (which may be one term or a sum of terms),
    antisymmetrizes over all factors that contain a down index (i.e. in a factor
    of the form something_index where ‘index’ is one or more letters).
    
    For example, given:
      "U1_i*P2^i*U2_c*U3^c*U1_j*P2^j*U2_d*U3^d"
    this function extracts the down-indexed factors U1_i, U2_c, U1_j, U2_d (in order),
    then forms the full antisymmetrization (i.e. the signed sum over all permutations of
    the index list [i, c, j, d]). Each permutation produces a new term in which the
    corresponding indices in the factors are replaced by the permuted ones.
    
    The new terms (with appropriate plus/minus signs) are joined with " + " and returned
    as a string.
    
    (Note: This function assumes that any factor with a down index is exactly of the form
    "[alphanumerics]_[letters]". It does not attempt to antisymmetrize up‐indices.)
    """
    # We'll split the equation on top-level '+' signs.
    # (Here, for simplicity, we assume the input is a single product term;
    #  if multiple terms are present, you may loop over them.)
    terms = [t.strip() for t in equation.split('+') if t.strip()]
    
    # A regex to recognize a factor with a down index.
    # For example: "U1_i" will yield base="U1" and index="i"
    down_index_pattern = re.compile(r'^([A-Za-z0-9]+)_([A-Za-z]+)$')
    
    def permutation_sign(original, permutation):
        """
        Compute the sign (+1 or -1) of the permutation that sends original --> permutation.
        We use a simple inversion count.
        """
        # Create a mapping from element to its position in the original order.
        pos_map = {v: i for i, v in enumerate(original)}
        # Map the permutation to the indices in the original ordering.
        perm_indices = [pos_map[x] for x in permutation]
        inv_count = 0
        n = len(perm_indices)
        for i in range(n):
            for j in range(i+1, n):
                if perm_indices[i] > perm_indices[j]:
                    inv_count += 1
        return -1 if inv_count % 2 else 1

    antisym_terms = []
    
    for term in terms:
        # Split term by '*' (assuming factors are separated by '*')
        factors = [f.strip() for f in term.split('*') if f.strip()]
        # Locate the positions (and record base and current index) of factors that have a down-index.
        indexed_positions = []
        for i, fac in enumerate(factors):
            m = down_index_pattern.match(fac)
            if m:
                base = m.group(1)
                idx = m.group(2)
                indexed_positions.append((i, base, idx))
        # If no down-indexed factors appear, keep the term unchanged.
        if not indexed_positions:
            antisym_terms.append(term)
        else:
            # Extract the original index list (in order of appearance).
            orig_idxs = [tup[2] for tup in indexed_positions]
            # Generate all distinct permutations of these indices.
            unique_perms = set(itertools.permutations(orig_idxs))
            # For each permutation, compute its sign and produce a new term.
            permuted_terms = []
            for perm in unique_perms:
                sgn = permutation_sign(orig_idxs, list(perm))
                # Build a new list of factors by replacing each down-indexed factor's index.
                new_factors = factors.copy()
                for (pos, base, old_idx), new_idx in zip(indexed_positions, perm):
                    new_factors[pos] = f"{base}_{new_idx}"
                # Form the product term (joined by '*')
                new_product = '*'.join(new_factors)
                # Prepend a '-' if sgn is -1.
                if sgn == -1:
                    new_product = '-' + new_product
                permuted_terms.append(new_product)
            # Join the permuted terms with " + " (representing the sum).
            antisym_terms.append(" + ".join(permuted_terms))
    
    # Finally, join all terms (if there were more than one top-level term) with " + "
    return " + ".join(antisym_terms)

def perform_full_operation(equation):
    while 'a1' in equation or 'a2' in equation or 'a3' in equation:
        equation = move_all_a_left(equation)
    equation = remove_traces(equation)
    equation = pull_all_non_canon_UP_left(equation)
    equation = move_all_a_left(equation)
    equation = pull_all_Divs_right(equation)
    equation = move_all_a_left(equation)
    equation = pull_all_PiPjs_left(equation)
    equation = move_all_a_left(equation)
    equation = pull_all_Bs_right(equation)
    equation = move_all_a_left(equation)
    equation = remove_traces(equation)
    return equation

def impose_y_z_rules(equation):
    equation = reorder_un_terms(equation)
    return impose_y_rule(impose_z_rule(equation))


gauge_variation_results = []
gauge_variation_z_results = []

for eq in rearrange_y_comms_all:
    while ('U' or 'P') in impose_y_z_rules(eq):
        eq = perform_full_operation(eq)
    gauge_variation_results.append(impose_y_z_rules(reorder_un_terms(eq)))

for eq in rearrange_y_comms_all_z:
    while ('U' or 'P') in impose_y_z_rules(eq):
        eq = perform_full_operation(eq)
    gauge_variation_z_results.append(impose_y_z_rules(reorder_un_terms(eq)))

full_gauge_variation_equation = ' + '.join(gauge_variation_results) + ' + ' + ' + '.join(gauge_variation_z_results)


two_deriv_DDIs_combs = [('P1_a', 'U1_b', 'U2_c', 'U3_e', 'P2^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'U1_b', 'U2_c', 'U3_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'U1_b', 'U2_c', 'U3_e', 'P1^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'U1_b', 'U2_c', 'U3_e', 'P2^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'U1_b', 'U2_c', 'U3_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'U1_b', 'U2_c', 'U3_e', 'P1^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'U1_b', 'U2_c', 'U3_e', 'P2^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'U1_b', 'U2_c', 'U3_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P2_a', 'U1_b', 'U2_c', 'U3_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P2_a', 'U1_b', 'U2_c', 'U3_e', 'P1^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P2_a', 'U1_b', 'U2_c', 'U3_e', 'P2^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P2_a', 'U1_b', 'U2_c', 'U3_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P2_a', 'U1_b', 'U2_c', 'U3_e', 'P1^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P2_a', 'U1_b', 'U2_c', 'U3_e', 'P2^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P2_a', 'U1_b', 'U2_c', 'U3_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P3_a', 'U1_b', 'U2_c', 'U3_e', 'P1^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P3_a', 'U1_b', 'U2_c', 'U3_e', 'P2^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P3_a', 'U1_b', 'U2_c', 'U3_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P3_a', 'U1_b', 'U2_c', 'U3_e', 'P1^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P3_a', 'U1_b', 'U2_c', 'U3_e', 'P2^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P3_a', 'U1_b', 'U2_c', 'U3_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'U1_b', 'U2_c', 'U3_e', 'P2^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'U1_b', 'U2_c', 'U3_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'U1_b', 'U2_c', 'U3_e', 'P1^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'U1_b', 'U2_c', 'U3_e', 'P2^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'U1_b', 'U2_c', 'U3_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P2_a', 'U1_b', 'U2_c', 'U3_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P2_a', 'U1_b', 'U2_c', 'U3_e', 'P1^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P2_a', 'U1_b', 'U2_c', 'U3_e', 'P2^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P2_a', 'U1_b', 'U2_c', 'U3_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P3_a', 'U1_b', 'U2_c', 'U3_e', 'P1^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P3_a', 'U1_b', 'U2_c', 'U3_e', 'P2^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P3_a', 'U1_b', 'U2_c', 'U3_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'U1_b', 'U2_c', 'U3_e', 'P2^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'U1_b', 'U2_c', 'U3_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P2_a', 'U1_b', 'U2_c', 'U3_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e')]
three_deriv_DDIs_combs = [('P1_a', 'P2_b', 'U1_c', 'U2_e', 'P1^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'P2_b', 'U1_c', 'U2_e', 'P2^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'P2_b', 'U1_c', 'U2_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'P2_b', 'U1_c', 'U3_e', 'P1^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'P2_b', 'U1_c', 'U3_e', 'P2^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'P2_b', 'U1_c', 'U3_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'P2_b', 'U2_c', 'U3_e', 'P1^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'P2_b', 'U2_c', 'U3_e', 'P2^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'P2_b', 'U2_c', 'U3_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'P1_b', 'U1_c', 'U2_e', 'P1^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'P1_b', 'U1_c', 'U2_e', 'P2^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'P1_b', 'U1_c', 'U2_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'P1_b', 'U1_c', 'U3_e', 'P1^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'P1_b', 'U1_c', 'U3_e', 'P2^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'P1_b', 'U1_c', 'U3_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'P1_b', 'U2_c', 'U3_e', 'P1^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'P1_b', 'U2_c', 'U3_e', 'P2^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'P1_b', 'U2_c', 'U3_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P1_a', 'U1_b', 'U2_c', 'U3_e', 'P2^a', 'P1^b', 'U1^c', 'U2^e'),
 ('P1_a', 'U1_b', 'U2_c', 'U3_e', 'P2^a', 'P1^b', 'U1^c', 'U3^e'),
 ('P1_a', 'U1_b', 'U2_c', 'U3_e', 'P2^a', 'P1^b', 'U2^c', 'U3^e'),
 ('P2_a', 'P1_b', 'U1_c', 'U2_e', 'P2^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P2_a', 'P1_b', 'U1_c', 'U2_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P2_a', 'P1_b', 'U1_c', 'U3_e', 'P2^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P2_a', 'P1_b', 'U1_c', 'U3_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P2_a', 'P1_b', 'U2_c', 'U3_e', 'P2^a', 'U1^b', 'U2^c', 'U3^e'),
 ('P2_a', 'P1_b', 'U2_c', 'U3_e', 'P3^a', 'U1^b', 'U2^c', 'U3^e')]
four_deriv_DDIs_combs = [('P1_a', 'P2_b', 'U2_c', 'U1_e', 'P1^a', 'P2^b', 'U2^c', 'U3^e'),
                         ('P2_a', 'P3_b', 'U3_c', 'U2_e', 'P2^a', 'P3^b', 'U3^c', 'U1^e'),
                         ('P3_a', 'P1_b', 'U1_c', 'U3_e', 'P3^a', 'P1^b', 'U1^c', 'U2^e'),
                         ('P3_a', 'P1_b', 'U1_c', 'U3_e', 'P1^a', 'P2^b', 'U3^c', 'U1^e'),
                         ('P1_a', 'P2_b', 'U2_c', 'U1_e', 'P2^a', 'P3^b', 'U1^c', 'U2^e'),
                         ('P2_a', 'P3_b', 'U3_c', 'U2_e', 'P3^a', 'P1^b', 'U2^c', 'U3^e')]


two_deriv_DDIs = [antisymmetrise_down_indices('*'.join(term)) for term in two_deriv_DDIs_combs]
three_deriv_DDIs = [antisymmetrise_down_indices('*'.join(term)) for term in three_deriv_DDIs_combs]
four_deriv_DDIs = [antisymmetrise_down_indices('*'.join(term)) for term in four_deriv_DDIs_combs]

naive_three_deriv_DDIs = [combine_same_terms(combine_powers_in_terms(fully_expand_equation(impose_dimension(impose_y_z_rules(perform_full_operation(perform_full_operation(perform_full_operation(perform_full_operation(equation))))))))) for equation in three_deriv_DDIs]
naive_two_deriv_DDIs = [combine_same_terms(combine_powers_in_terms(fully_expand_equation(impose_dimension(impose_y_z_rules(perform_full_operation(perform_full_operation(perform_full_operation(perform_full_operation(equation))))))))) for equation in two_deriv_DDIs]
naive_four_deriv_DDIs = [combine_same_terms(combine_powers_in_terms(fully_expand_equation(impose_dimension(impose_y_z_rules(perform_full_operation(perform_full_operation(perform_full_operation(perform_full_operation(equation))))))))) for equation in four_deriv_DDIs]

naive_three_deriv_DDIs_flat = [[term for term in top_level_split(equation) if 'l^2' not in term and 'm' not in term] for equation in naive_three_deriv_DDIs]
naive_two_deriv_DDIs_flat = [[term for term in top_level_split(equation) if 'l^2' not in term and 'm' not in term] for equation in naive_two_deriv_DDIs]
naive_four_deriv_DDIs_flat = [[term for term in top_level_split(equation) if 'l^2' not in term and 'm' not in term] for equation in naive_four_deriv_DDIs]

G_minus_Y1Z1_sqr_DDI = flip_eq_sign(two_deriv_DDIs[5])
G_minus_Y2Z2_sqr_DDI = flip_eq_sign(two_deriv_DDIs[10])
G_minus_Y2Z2_sqr_DDI = flip_eq_sign(two_deriv_DDIs[17])

Y1Z1G_DDI = two_deriv_DDIs[11]
Y2Z2G_DDI = two_deriv_DDIs[1]
Y3Z3G_DDI = two_deriv_DDIs[0]

Y1Y3sqrZ3_DDI = three_deriv_DDIs[3]
Y2Y3sqrZ3_DDI = flip_eq_sign(three_deriv_DDIs[7])
Y1Y2sqrZ2_DDI = flip_eq_sign(three_deriv_DDIs[0])
Y3Y2sqrZ2_DDI = flip_eq_sign(three_deriv_DDIs[8])
Y2Y1sqrZ1_DDI = flip_eq_sign(three_deriv_DDIs[1])
Y3Y1sqrZ1_DDI = three_deriv_DDIs[5]

Y1sqrY2Y3_DDI = four_deriv_DDIs[2]
Y1Y2sqrY3_DDI = four_deriv_DDIs[0]
Y1Y2Y3sqr_DDI = four_deriv_DDIs[1]

# general DDIs


# two derivatives
DDI_1 = '(' + two_deriv_DDIs[11] + ')*Z1^n1*Z2^n2*Z3^n3'
DDI_2 = '(' + two_deriv_DDIs[1] + ')*Z1^n1*Z2^n2*Z3^n3'       # antisym*V(z)
DDI_3 = '(' + two_deriv_DDIs[0] + ')*Z1^n1*Z2^n2*Z3^n3'
DDI_4 = '(' + flip_eq_sign(two_deriv_DDIs[5]) + ')*Z1^n1*Z2^n2*Z3^n3'
DDI_5 = '(' + flip_eq_sign(two_deriv_DDIs[10]) + ')*Z1^n1*Z2^n2*Z3^n3'       # antisym*V(z)
DDI_6 = '(' + flip_eq_sign(two_deriv_DDIs[17]) + ')*Z1^n1*Z2^n2*Z3^n3'

# three derivatives
DDI_7 = '(' + three_deriv_DDIs[3] + ')*Z1^n1*Z2^n2*Z3^n3'
DDI_8 = '(' + flip_eq_sign(three_deriv_DDIs[7]) + ')*Z1^n1*Z2^n2*Z3^n3'
DDI_9 = '(' + flip_eq_sign(three_deriv_DDIs[0]) + ')*Z1^n1*Z2^n2*Z3^n3'
DDI_10 = '(' + flip_eq_sign(three_deriv_DDIs[8]) + ')*Z1^n1*Z2^n2*Z3^n3'      # antisym*V(z)
DDI_11 = '(' + flip_eq_sign(three_deriv_DDIs[1]) + ')*Z1^n1*Z2^n2*Z3^n3'
DDI_12 = '(' + three_deriv_DDIs[5] + ')*Z1^n1*Z2^n2*Z3^n3'
DDI_13 = '(' + flip_eq_sign(two_deriv_DDIs[5]) + ')*Y3*Z1^n1*Z2^n2*Z3^n3'
DDI_14 = '(' + flip_eq_sign(two_deriv_DDIs[10]) + ')*Y1*Z1^n1*Z2^n2*Z3^n3'
DDI_15 = '(' + flip_eq_sign(two_deriv_DDIs[17]) + ')*Y2*Z1^n1*Z2^n2*Z3^n3'

# four derivatives
DDI_16 = '(' + four_deriv_DDIs[2] + ')*Z1^n1*Z2^n2*Z3^n3'
DDI_17 = '(' + four_deriv_DDIs[0] + ')*Z1^n1*Z2^n2*Z3^n3'
DDI_18 = '(' + four_deriv_DDIs[1] + ')*Z1^n1*Z2^n2*Z3^n3'
DDI_19 = '(' + four_deriv_DDIs[4] + ')*Z1^n1*Z2^n2*Z3^n3'
DDI_20 = '(' + four_deriv_DDIs[5] + ')*Z1^n1*Z2^n2*Z3^n3'
DDI_21 = '(' + four_deriv_DDIs[3] + ')*Z1^n1*Z2^n2*Z3^n3'
DDI_22 = '(' + three_deriv_DDIs[3] + ')*Y3*Z1^n1*Z2^n2*Z3^n3'
DDI_23 = '(' + flip_eq_sign(three_deriv_DDIs[7]) + ')*Y3*Z1^n1*Z2^n2*Z3^n3'         # DDIs cyclic 7 -> 0 -> 5
DDI_24 = '(' + flip_eq_sign(three_deriv_DDIs[0]) + ')*Y2*Z1^n1*Z2^n2*Z3^n3'         # and 3 -> 8 -> 1
DDI_25 = '(' + flip_eq_sign(three_deriv_DDIs[8]) + ')*Y2*Z1^n1*Z2^n2*Z3^n3'
DDI_26 = '(' + flip_eq_sign(three_deriv_DDIs[1]) + ')*Y1*Z1^n1*Z2^n2*Z3^n3'
DDI_27 = '(' + three_deriv_DDIs[5] + ')*Y1*Z1^n1*Z2^n2*Z3^n3'
DDI_28 = '(' + flip_eq_sign(two_deriv_DDIs[5]) + ')*Y2*Y2*Z1^n1*Z2^n2*Z3^n3'
DDI_29 = '(' + flip_eq_sign(two_deriv_DDIs[10]) + ')*Y3*Y3*Z1^n1*Z2^n2*Z3^n3'
DDI_30 = '(' + flip_eq_sign(two_deriv_DDIs[17]) + ')*Y1*Y1*Z1^n1*Z2^n2*Z3^n3'

all_DDIs = [DDI_1, DDI_2, DDI_3, DDI_4, DDI_5, DDI_6, DDI_7, DDI_8, DDI_9, DDI_10,
             DDI_11, DDI_12, DDI_13, DDI_14, DDI_15, DDI_16, DDI_17, DDI_18, DDI_19, 
             DDI_20, DDI_21, DDI_22, DDI_23, DDI_24, DDI_25, DDI_26, DDI_27, DDI_28, 
             DDI_29, DDI_30]

all_DDIs_processed = []
for DDI in all_DDIs:
    eq = replace_all_Ys(DDI)
    eq = fully_expand_equation(eq)
    while ('U' or 'P') in impose_y_z_rules(eq):
        eq = perform_full_operation(eq)
        print(all_DDIs.index(DDI))
    eq = impose_y_z_rules(eq)
    eq = impose_dimension(eq)
    eq = fully_expand_equation(eq)
    eq = combine_powers_in_terms(eq)
    eq = combine_same_terms(eq)
    all_DDIs_processed.append(eq)

def reduce_linear_dependencies_DDIs(DDI_list):
    for DDI_1 in DDI_list:
        highest_deriv_terms_in_DDI_1 = [term for term in top_level_split(DDI_1) if 'l' not in term and 'm' not in term]
        DDI_1_pos = DDI_list.index(DDI_1)
        for DDI_2 in DDI_list[DDI_1_pos+1:]:
            highest_deriv_terms_in_DDI_2 = [term for term in top_level_split(DDI_2) if 'l' not in term and 'm' not in term]
            DDI_2_pos = DDI_list.index(DDI_2)
            if sorted(highest_deriv_terms_in_DDI_1) == sorted(highest_deriv_terms_in_DDI_2) and DDI_1_pos != DDI_2_pos:
                DDI_list[DDI_1_pos] = combine_same_terms(DDI_1 + ' + ' + flip_eq_sign(DDI_2))

    for DDI in DDI_list:
        all_terms_have_l = True
        for term in top_level_split(DDI):
            if 'l' not in term:
                all_terms_have_l = False
                break
        if all_terms_have_l:
            terms_removed_l = []
            for term in top_level_split(DDI):
                terms_removed_l.append(term.replace('(1/l^2)', '', 1))
            DDI_removed_l = '*'.join(terms_removed_l)
            DDI_list[DDI_list.index(DDI)] = DDI_removed_l
    return DDI_list

def rotate_expr(expr, shift):
    """
    Rotates the indices in the expression by a given shift.
    For example, with shift=1:
      n1 -> n2, n2 -> n3, n3 -> n1,
      Z1 -> Z2, Z2 -> Z3, Z3 -> Z1,
      Y1 -> Y2, Y2 -> Y3, Y3 -> Y1,
      s1 -> s2, s2 -> s3, s3 -> s1.
    """
    # Build a mapping for the digits under a cyclic shift.
    mapping = {str(i): str(((i - 1 + shift) % 3) + 1) for i in range(1, 4)}
    # This regex finds occurrences of n, Z, Y, or s followed by a digit.
    pattern = re.compile(r'([nZYsmUP])(\d)')
    
    def repl(match):
        var = match.group(1)
        digit = match.group(2)
        return var + mapping[digit]
    
    return pattern.sub(repl, expr)

def canonicalize_term(term):
    """
    Canonicalizes a single term by standardizing the order of factors.
    Assumes factors are separated by '*' (multiplication).
    Preserves a leading minus sign if present.
    """
    term = term.strip()
    sign = ""
    if term.startswith("-"):
        sign = "-"
        term = term[1:].strip()
    # Split the term by '*' and remove empty factors.
    factors = [f.strip() for f in term.split("*") if f.strip()]
    # Sort factors lexicographically.
    factors.sort()
    return sign + "*".join(factors)

def canonicalize_terms(expr):
    """
    Canonicalizes an entire expression by:
      - Normalizing whitespace.
      - Splitting into summation terms (assuming terms are delimited by " + ").
      - Canonicalizing each term (ordering the factors).
      - Sorting the canonicalized terms.
    """
    expr = " ".join(expr.split())
    # Split on " + " (assumes that each term is separated by " + ").
    terms = expr.split(" + ")
    canonical_terms = [canonicalize_term(term) for term in terms if term.strip()]
    canonical_terms.sort()
    return " + ".join(canonical_terms)

def canonical(expr):
    """
    Returns a canonical representation for an expression by considering all
    cyclic rotations (0, 1, and 2) and applying term and factor order standardization.
    """
    # Generate all rotated forms.
    rotations = [rotate_expr(expr, shift) for shift in range(3)]
    # For each rotated form, canonicalize the order of summation terms and factors.
    rotations = [canonicalize_terms(rot) for rot in rotations]
    # Return the lexicographically smallest representation.
    return min(rotations)

def filter_expressions(expr_list):
    """
    Given a list of expressions, returns a filtered list containing only one representative
    for each group of expressions that are equivalent up to cyclic permutation of indices,
    ordering of summation terms, and ordering of factors within each term.
    """
    unique = {}
    for expr in expr_list:
        key = canonical(expr)
        if key not in unique:
            unique[key] = expr  # or store key if you want the canonical version
    return list(unique.values())

filtered_list = filter_expressions(all_DDIs_processed)

def extract_relevant_Z_power(equation, power=1):
    terms = top_level_split(equation)
    new_terms = []
    for term in terms:
        sub_terms_list = term.split('*')
        if 'Y3^2' in term:
            if '-Y3^2' in term:
                Y3sqr_pos = sub_terms_list.index('-Y3^2')
            else:
                Y3sqr_pos = sub_terms_list.index('Y3^2')
            Z3_var = [var for var in sub_terms_list if var.startswith('Z3')][0]
            Z3_pos = sub_terms_list.index(Z3_var)
            Z3_exponent = Z3_var.split('^')[1]
            for i in range(power):
                Z3_exponent = decrement_exponent(Z3_exponent)
            sub_terms_list[Z3_pos] = 'Z3^' + Z3_exponent
            if power == 1:
                sub_terms_list.insert(Y3sqr_pos+1, 'Z3')
            else:
                sub_terms_list.insert(Y3sqr_pos+1, f'Z3^{power}')
        if 'Y2^2' in term:
            if '-Y2^2' in term:
                Y2sqr_pos = sub_terms_list.index('-Y2^2')
            else:
                Y2sqr_pos = sub_terms_list.index('Y2^2')
            Z2_var = [var for var in sub_terms_list if var.startswith('Z2')][0]
            Z2_pos = sub_terms_list.index(Z2_var)
            Z2_exponent = Z2_var.split('^')[1]
            for i in range(power):
                Z2_exponent = decrement_exponent(Z2_exponent)
            sub_terms_list[Z2_pos] = 'Z2^' + Z2_exponent
            if power == 1:
                sub_terms_list.insert(Y2sqr_pos+1, 'Z2')
            else:
                sub_terms_list.insert(Y2sqr_pos+1, f'Z2^{power}')
        if 'Y1^2' in term:
            if '-Y1^2' in term:
                Y1sqr_pos = sub_terms_list.index('-Y1^2')
            else:
                Y1sqr_pos = sub_terms_list.index('Y1^2')
            Z1_var = [var for var in sub_terms_list if var.startswith('Z1')][0]
            Z1_pos = sub_terms_list.index(Z1_var)
            Z1_exponent = Z1_var.split('^')[1]
            for i in range(power):
                Z1_exponent = decrement_exponent(Z1_exponent)
            sub_terms_list[Z1_pos] = 'Z1^' + Z1_exponent
            if power == 1:
                sub_terms_list.insert(Y1sqr_pos+1, 'Z1')
            else:
                sub_terms_list.insert(Y1sqr_pos+1, f'Z1^{power}')
        new_terms.append('*'.join(sub_terms_list))
    
    return ' + '.join(reversed(new_terms))

def DDI_sub_Y_cubed(equation):
    # want to find terms with Yi^3 and substitute to bring into Yi^2*Yj or Y1*Y2*Y3 form
    equation = re.sub(r'Y3\^3', 'Y3*Y3^2', equation)
    equation = re.sub(r'Y2\^3', 'Y2*Y2^2', equation)
    equation = re.sub(r'Y1\^3', 'Y1*Y1^2', equation)

    equation = extract_relevant_Z_power(equation, power=2)

    anti_sym1 = '(' + 'U1_a*P2^a*U1_b*P2^b*U2_c*U3^c*U2_d*U3^d' + ' + ' + Y1Z1G_DDI + ')'
    anti_sym2 = '(' + 'U2_a*P3^a*U2_b*P3^b*U3_c*U1^c*U3_d*U1^d' + ' + ' + Y2Z2G_DDI + ')'
    anti_sym3 = '(' + 'U3_a*P1^a*U3_b*P1^b*U1_c*U2^c*U1_d*U2^d' + ' + ' + Y3Z3G_DDI + ')'

    equation = re.sub(r'Y1\^2\*Z1\^2', anti_sym1, equation)
    equation = re.sub(r'Y2\^2\*Z2\^2', anti_sym2, equation)
    equation = re.sub(r'Y3\^2\*Z3\^2', anti_sym3, equation)

    terms = top_level_split(equation)
    new_terms = []
    for term in terms:
        term = replace_Y1_with_U1P2(replace_Y2_with_U2P3(replace_Y3_with_U3P1(term)))
        new_terms.append(term)

    return combine_powers_in_terms(fully_expand_equation(' + '.join(new_terms)))

def DDI_sub_YZG(equation):
    expanded_Y1Z1G = replace_all_YZs(fully_expand_equation('Y1*Z1*(Y1*Z1+Y2*Z2+Y3*Z3)'))
    expanded_Y2Z2G = replace_all_YZs(fully_expand_equation('Y2*Z2*(Y1*Z1+Y2*Z2+Y3*Z3)'))
    expanded_Y3Z3G = replace_all_YZs(fully_expand_equation('Y3*Z3*(Y1*Z1+Y2*Z2+Y3*Z3)'))

    anti_sym1 = '(' + expanded_Y1Z1G + ' + ' + Y1Z1G_DDI + ')'
    anti_sym2 = '(' + expanded_Y2Z2G + ' + ' + Y2Z2G_DDI + ')'
    anti_sym3 = '(' + expanded_Y3Z3G + ' + ' + Y3Z3G_DDI + ')'

    equation = re.sub(r'Y1\*Z1\*\(Y1\*Z1\+Y2\*Z2\+Y3\*Z3\)', anti_sym1, equation)
    equation = re.sub(r'Y2\*Z2\*\(Y1\*Z1\+Y2\*Z2\+Y3\*Z3\)', anti_sym2, equation)
    equation = re.sub(r'Y3\*Z3\*\(Y1\*Z1\+Y2\*Z2\+Y3\*Z3\)', anti_sym3, equation)

    return pull_minus_signs_to_front(equation)

def DDI_sub_YYYZ(equation):
    Y1Y3sqr = '(' + 'U1_w*P2^w*U3_x*P1^x*U3_y*P1^y*U1_z*U2^z' + ' + ' + Y1Y3sqrZ3_DDI + ')'
    Y2Y3sqr = '(' + 'U2_w*P3^w*U3_x*P1^x*U3_y*P1^y*U1_z*U2^z' + ' + ' + Y2Y3sqrZ3_DDI + ')'
    Y1Y2sqr = '(' + 'U1_w*P2^w*U2_x*P3^x*U2_y*P3^y*U1_z*U3^z' + ' + ' + Y1Y2sqrZ2_DDI + ')'
    Y3Y2sqr = '(' + 'U3_w*P1^w*U2_x*P3^x*U2_y*P3^y*U1_z*U3^z' + ' + ' + Y3Y2sqrZ2_DDI + ')'
    Y2Y1sqr = '(' + 'U2_w*P3^w*U1_x*P2^x*U1_y*P2^y*U2_z*U3^z' + ' + ' + Y2Y1sqrZ1_DDI + ')'
    Y3Y1sqr = '(' + 'U3_w*P1^w*U1_x*P2^x*U1_y*P2^y*U2_z*U3^z' + ' + ' + Y3Y1sqrZ1_DDI + ')'

    equation = extract_relevant_Z_power(equation)

    equation = re.sub(r'\bY1\*Y3\^2\*Z3\b', Y1Y3sqr, equation)
    equation = re.sub(r'\bY2\*Y3\^2\*Z3\b', Y2Y3sqr, equation)
    equation = re.sub(r'\bY1\*Y2\^2\*Z2\b', Y1Y2sqr, equation)
    equation = re.sub(r'\bY2\^2\*Z2\*Y3\b', Y3Y2sqr, equation)
    equation = re.sub(r'\bY1\^2\*Z1\*Y2\b', Y2Y1sqr, equation)
    equation = re.sub(r'\bY1\^2\*Z1\*Y3\b', Y3Y1sqr, equation)

    return fully_expand_equation(equation)

def DDI_sub_YYYY(equation):
    Y1sqrY2Y3 = '(' + 'U1_w*P2^w*U1_x*P2^x*U2_y*P3^y*U3_z*P1^z' + ' + ' + Y1sqrY2Y3_DDI + ')'
    Y1Y2sqrY3 = '(' + 'U2_w*P3^w*U2_x*P3^x*U3_y*P1^y*U1_z*P2^z' + ' + ' + Y1Y2sqrY3_DDI + ')'
    Y1Y2Y3sqr = '(' + 'U3_w*P1^w*U3_x*P1^x*U1_y*P2^y*U2_z*P3^z' + ' + ' + Y1Y2Y3sqr_DDI + ')'

    equation = re.sub(r'\bY1\^2\*Y2\*Y3\b', Y1sqrY2Y3, equation)
    equation = re.sub(r'\bY1\*Y2\^2\*Y3\b', Y1Y2sqrY3, equation)
    equation = re.sub(r'\bY1\*Y2\*Y3\^2\b', Y1Y2Y3sqr, equation)

    return fully_expand_equation(equation)

def fully_process_three_deriv_gauge_variation_eq(equation, no_deriv=2):
    equation = complete_sub_and_expansion(equation)
    if no_deriv == 2:
        equation = DDI_sub_Y_cubed(equation)
        equation = perform_full_operation(perform_full_operation(perform_full_operation(perform_full_operation(equation))))
        equation = impose_y_z_rules(equation)
        equation = complete_sub_and_expansion(equation)
        equation = DDI_sub_YYYZ(equation)
        equation = perform_full_operation(perform_full_operation(perform_full_operation(perform_full_operation(equation))))
        equation = impose_y_z_rules(equation)
        equation = complete_sub_and_expansion(equation)
    if no_deriv == 3:
        equation = DDI_sub_YYYY(equation)
        equation = perform_full_operation(perform_full_operation(perform_full_operation(perform_full_operation(equation))))
        equation = impose_y_z_rules(equation)
        equation = complete_sub_and_expansion(equation)
        equation = DDI_sub_Y_cubed(equation)
        equation = perform_full_operation(perform_full_operation(perform_full_operation(perform_full_operation(equation))))
        equation = impose_y_z_rules(equation)
        equation = complete_sub_and_expansion(equation)
        equation = DDI_sub_YYYZ(equation)
        equation = perform_full_operation(perform_full_operation(perform_full_operation(perform_full_operation(equation))))
        equation = impose_y_z_rules(equation)
        equation = complete_sub_and_expansion(equation)
    return equation

def extract_Y1_coeffs(equation):
    terms = top_level_split(equation)
    coeffs = []
    for term in terms:
        pattern = re.compile(r'Y1\*Z1\^\(n1\+2\)\*Z2\^\(n2\)\*Z3\^\(n3\)')
        if bool(pattern.search(term)):
            # Replace all occurrences with an empty string.
            coeff = pattern.sub("", term)
            coeffs.append(coeff.strip('*'))
        else: continue
    return ' + '.join(coeffs)

def extract_Y2_coeffs(equation):
    terms = top_level_split(equation)
    coeffs = []
    for term in terms:
        pattern = re.compile(r'Y2\*Z1\^\(n1\+1\)\*Z2\^\(n2\+1\)\*Z3\^\(n3\)')
        if bool(pattern.search(term)):
            # Replace all occurrences with an empty string.
            coeff = pattern.sub("", term)
            coeffs.append(coeff.strip('*'))
        else: continue
    return ' + '.join(coeffs)

def extract_Y3_coeffs(equation):

    terms = top_level_split(equation)
    coeffs = []
    for term in terms:
        pattern = re.compile(r'Y3\*Z1\^\(n1\+1\)\*Z2\^\(n2\)\*Z3\^\(n3\+1\)')
        if bool(pattern.search(term)):
            # Replace all occurrences with an empty string.
            coeff = pattern.sub("", term)
            coeffs.append(coeff.strip('*'))
        else: continue
    return ' + '.join(coeffs)

def extract_Y1Y3_coeffs(equation):

    terms = top_level_split(equation)
    coeffs = []
    for term in terms:
        pattern = re.compile(r'Y1\*Y3\*Z1\^\(p1\+1\)\*Z2\^\(p2\-1\)\*Z3\^\(p3\)')
        if bool(pattern.search(term)):
            # Replace all occurrences with an empty string.
            coeff = pattern.sub("", term)
            coeffs.append(coeff.strip('*'))
        else: continue
    return ' + '.join(coeffs)

def extract_Y1Y2_coeffs(equation):

    terms = top_level_split(equation)
    coeffs = []
    for term in terms:
        pattern = re.compile(r'Y1\*Y2\*Z1\^\(p1\+1\)\*Z2\^\(p2\)\*Z3\^\(p3\-1\)')
        if bool(pattern.search(term)):
            # Replace all occurrences with an empty string.
            coeff = pattern.sub("", term)
            coeffs.append(coeff.strip('*'))
        else: continue
    return ' + '.join(coeffs)

def extract_Y2Y3_coeffs(equation):

    terms = top_level_split(equation)
    coeffs = []
    for term in terms:
        pattern = re.compile(r'Y2\*Y3\*Z1\^\(p1\)\*Z2\^\(p2\)\*Z3\^\(p3\)')
        if bool(pattern.search(term)):
            # Replace all occurrences with an empty string.
            coeff = pattern.sub("", term)
            coeffs.append(coeff.strip('*'))
        else: continue
    return ' + '.join(coeffs)

def extract_Y1_sqr_coeffs(equation):

    terms = top_level_split(equation)
    coeffs = []
    for term in terms:
        pattern = re.compile(r'Y1\^2\*Z1\^\(p1\+1\)\*Z2\^\(p2\-1\)\*Z3\^\(p3\)')
        if bool(pattern.search(term)):
            # Replace all occurrences with an empty string.
            coeff = pattern.sub("", term)
            coeffs.append(coeff.strip('*'))
        else: continue
    return ' + '.join(coeffs)

def extract_Y2_sqr_coeffs(equation):

    terms = top_level_split(equation)
    coeffs = []
    for term in terms:
        pattern = re.compile(r'Y2\^2\*Z1\^\(p1\)\*Z2\^\(p2\+1\)\*Z3\^\(p3\-1\)')
        if bool(pattern.search(term)):
            # Replace all occurrences with an empty string.
            coeff = pattern.sub("", term)
            coeffs.append(coeff.strip('*'))
        else: continue
    return ' + '.join(coeffs)

def extract_Y3_sqr_coeffs(equation):

    terms = top_level_split(equation)
    coeffs = []
    for term in terms:
        pattern = re.compile(r'Y3\^2\*Z1\^\(p1\)\*Z2\^\(p2\-1\)\*Z3\^\(p3\+1\)')
        if bool(pattern.search(term)):
            # Replace all occurrences with an empty string.
            coeff = pattern.sub("", term)
            coeffs.append(coeff.strip('*'))
        else: continue
    return ' + '.join(coeffs)

def extract_Z1Z2Z3_coeffs(equation):

    terms = top_level_split(equation)
    coeffs = []
    for term in terms:
        pattern = re.compile(r'Z1\^\(p1\+1\)\*Z2\^\(p2\)\*Z3\^\(p3\)')
        if bool(pattern.search(term)):
            # Replace all occurrences with an empty string.
            coeff = pattern.sub("", term)
            coeffs.append(coeff.strip('*'))
        else: continue
    return ' + '.join(coeffs)

final_gauge_variation_equation = fully_process_three_deriv_gauge_variation_eq(full_gauge_variation_equation, no_deriv=2)

A_red_sym = '(1/2)*(1/l^2)*(2 + s + s + 2 - 2*s^2 - 2*s*2 - 2*s*2 + s^3 + 4*s + s^3 + 4*s + 2*s^2 + 2*s^2 - 4*s^2 - s^3 - s^3 - 8)'

A_sym = '(1/2)*(1/l^2)*(2 + s1 + s2 + s3 - 2*s1*s2 - 2*s1*s3 - 2*s2*s3 + s1*s2^2 + s1*s3^2 + s2*s1^2 + s2*s3^2 + s3*s1^2 + s3*s2^2 - 2*s1*s2*s3 - s1^3 - s2^3 - s3^3)'

A_red_non_sym = '(1/2)*(1/l^2)*(6 + s - 3*s - 7*2 - 4*s^2 - 4*s^2 + 4*4 + 6*s^2 + 4*s*2 + s^3 + s^3 - 8 - s^3 - s*4 - s^3 - s*4 + 2*s^2 + 2*s^2 - 2*s^2*2)'

A_non_sym = '(1/2)*(1/l^2)*(6 + s1 - 3*s2 - 7*s3 - 4*s1^2 - 4*s2^2 + 4*s3^2 + 6*s1*s2 + 4*s2*s3 + s1^3 + s2^3 - s3^3 - s1*s2^2 - s1*s3^2 - s2*s1^2 - s2*s3^2 + s3*s1^2 + s3*s2^2 - 2*s1*s2*s3)'

# three derivative vertex gauge variation constraints:

C_cond = '(B + (1/2)*(1/l^2)*(3*s2 + -s2^2 + -5*s3 + s3^2 + 2*s1 + -s1*s2 + s1*s3))'

A_cond = combine_same_terms(extract_Y1Y3_coeffs(final_gauge_variation_equation) + ' + ' + extract_Y1Y2_coeffs(final_gauge_variation_equation))

A = '(1/2)*(1/l^2)*(2 + -2*s3 - 5*s1 + s1*s2 + s1*s3 + s1^2)'
B = '(1/2)*(1/l^2)*(2 + -2*s1 - 5*s2 + s2*s3 + s2*s1 + s2^2)'
C = '(1/2)*(1/l^2)*(2 + -2*s2 - 5*s3 + s3*s1 + s3*s2 + s3^2)'


end_time = time.perf_counter()
print(f"Elapsed time: {end_time - start_time} seconds")