import re
import time
from collections import Counter
from collections import defaultdict
from Gauge_Variation import filtered_list
from fractions import Fraction
import itertools


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

def write_Y_powers_out(equation):
    terms = top_level_split(equation)
    resulting_terms = []
    for term in terms:
        sub_terms_list = term.split('*')
        for var in sub_terms_list:
            if (var.startswith('Y') or var.startswith('-Y')) and '^' in var:
                Y_var, Y_power = var.split('^')
                if Y_var[0] == '-':
                    Y_var = Y_var[1:]
                    minus = True
                else:
                    minus = False
                new_Y_string = []
                for i in range(int(Y_power)):
                    new_Y_string.append(Y_var)
                if minus:
                    new_Y_string = '-' + '*'.join(new_Y_string)
                else:
                    new_Y_string = '*'.join(new_Y_string)
                sub_terms_list[sub_terms_list.index(var)] = new_Y_string
        resulting_terms.append('*'.join(sub_terms_list))
    return ' + '.join(resulting_terms)

def flip_eq_sign(equation):
    terms = top_level_split(equation)
    sign_flipped_terms = []
    for term in terms:
        if term[0] == '-':
            sign_flipped_terms.append(term[1:])
        else:
            sign_flipped_terms.append('-' + term)
        
    return ' + '.join(sign_flipped_terms)

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

all_DDIs_processed = filtered_list

def reduce_linear_dependencies_DDIs(DDI_list):
    for i, DDI_1 in enumerate(DDI_list):
        highest_deriv_terms_in_DDI_1 = [term for term in top_level_split(DDI_1) if 'l' not in term and 'm' not in term]
        for j in range(i + 1, len(DDI_list)):
            DDI_2 = DDI_list[j]
            highest_deriv_terms_in_DDI_2 = [term for term in top_level_split(DDI_2) if 'l' not in term and 'm' not in term]
            if sorted(highest_deriv_terms_in_DDI_1) == sorted(highest_deriv_terms_in_DDI_2) and i != j and len(highest_deriv_terms_in_DDI_1) != 0:
                combined = combine_same_terms(DDI_1 + ' + ' + flip_eq_sign(DDI_2))
                DDI_list[i] = combined
            elif sorted(highest_deriv_terms_in_DDI_1) == sorted(flip_eq_sign(highest_deriv_terms_in_DDI_2)):
                combined = combine_same_terms(DDI_1 + ' + ' + DDI_2)
                DDI_list[i] = combined

    for idx, DDI in enumerate(DDI_list):
        all_terms_have_l = True
        for term in top_level_split(DDI):
            if 'l' not in term:
                all_terms_have_l = False
                break
        if all_terms_have_l:
            terms_removed_l = []
            for term in top_level_split(DDI):
                # Remove one instance of '(1/l^2)' and tidy up
                cleaned = term.replace('(1/l^2)', '', 1).strip('*').replace('**', '*').replace('-*', '-')
                terms_removed_l.append(cleaned)
            DDI_removed_l = '*'.join(terms_removed_l)
            DDI_list[idx] = DDI_removed_l
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

'''
# DDI_40 should be linear dependent on other DDIs due to cyclic nature (53, 43, 40 generated from cyclic permutations)
filter_DDI_40 = filtered_list[filtered_list.index(all_DDIs_processed[39])]
filter_DDI_43 = filtered_list[filtered_list.index(all_DDIs_processed[42])]
# same for DDI 38 (51,48,38 generated from cyclic permutations)
filter_DDI_38 = filtered_list[filtered_list.index(all_DDIs_processed[37])]
filter_DDI_48 = filtered_list[filtered_list.index(all_DDIs_processed[47])]

# first subtract the 'cyclic' part e.g. some cyclic permutation of 53 or 43 for DDI_40 and 51 or 48 for 38 
filter_DDI_40 = combine_same_terms(filter_DDI_40 + ' + ' + flip_eq_sign(rotate_expr(filter_DDI_43,-1)))
filter_DDI_38 = combine_same_terms(filter_DDI_38 + ' + ' + flip_eq_sign(rotate_expr(filter_DDI_48,-1)))

# useful for finding which DDI to use for linear reduction
def obtain_filter_DDIs_with_n_derivs(n_derivs):
    return [[[term for term in top_level_split(eq) if 'l' not in term and 'm' not in term],filtered_list.index(eq)] for eq in filtered_list if len([term for term in top_level_split(eq) if 'l' not in term and 'm' not in term]) == n_derivs]

def get_all_terms_for_given_number_of_Ys(eq, n_Ys=2):
    terms_with_n_Ys = []
    for term in top_level_split(eq):
        number_of_Ys = len([1 for Y in write_Y_powers_out(term) if Y=='Y'])
        if number_of_Ys == n_Ys:
            terms_with_n_Ys.append(term)
    return terms_with_n_Ys

# Two derivative DDI for cancelling in 38
filter_DDI_21 = filtered_list[21]
# Two derivative DDI for cancelling in 40
filter_DDI_0 = filtered_list[0]

# DDI_38 is linearly dependent on the other DDIs under cyclic permutation -> can drop
filter_DDI_38 = combine_same_terms(filter_DDI_38 + ' + ' + fully_expand_equation('(1/l^2)*(' + rotate_expr(filter_DDI_21,1) + ')'))

# DDI_40 is linearly dependent on the other DDIs under cyclic permutation -> can drop
filter_DDI_40 = combine_same_terms(filter_DDI_40 + ' + ' + fully_expand_equation('-(1/l^2)*(' + rotate_expr(filter_DDI_0,0) + ')'))

# drop these
filtered_list[filtered_list.index(all_DDIs_processed[39])] = '0'
filtered_list[filtered_list.index(all_DDIs_processed[37])] = '0'
filtered_list = [eq for eq in filtered_list if eq != '0']
'''
'''
# now rearrange some of the DDIs - possibly remove or reduce lower derivative terms for 4 deriv DDIs
# first looks at entry 4

filter_DDI_4 = filtered_list[4]
# contains terms '(1/l^2)*(-7*Y1*Y2*Z2^2 + -4*Y1*Y3*Z2*Z3 + 2*Y1^2*Z1*Z2)*Z1^(n1)*Z2^(n2)*Z3^(n3+1)' and
# '2*n1*(1\l^2)*(-Y1*Y2*Z1*Z2 + -2*Y1*Y3*Z1*Z3 + -2*Y1^2*Z1^2)*Z1^(n1-1)*Z2^(n2+1)*Z3^(n3+1)'

filter_DDI_4 = combine_same_terms(filter_DDI_4 + ' + ' + fully_expand_equation('-2*n1*(1/l^2)*(' + rotate_expr(filter_DDI_0,0) + ')'))'
'''

# Want to automatically write DDIs in nicer factored/differential operators form

def pull_out_Z_factor(equation):
    terms = top_level_split(equation)
    new_terms = []
    for term in terms:
        variables = top_level_split(term, delimiter='*')
        new_vars = []
        for var in variables:
            if var.startswith('Z1^(n1+'):
                z_power = var.split('+')[1][:-1]
                if z_power == '1':
                    new_z1_term = 'Z1'
                else:
                    new_z1_term = f'Z1^{z_power}'
                new_vars.append(new_z1_term)
            elif var.startswith('Z1^(n1-'):
                z_power = var.split('-')[1][:-1]
                new_z1_term = f'Z1^(-{z_power})'
                new_vars.append(new_z1_term)
            elif var == 'Z1^(n1)' or var == 'Z1^n1':
                continue

            elif var.startswith('Z2^(n2+'):
                z_power = var.split('+')[1][:-1]
                if z_power == '1':
                    new_z2_term = 'Z2'
                else:
                    new_z2_term = f'Z2^{z_power}'
                new_vars.append(new_z2_term)
            elif var.startswith('Z2^(n2-'):
                z_power = var.split('-')[1][:-1]
                new_z2_term = f'Z2^(-{z_power})'
                new_vars.append(new_z2_term)
            elif var == 'Z2^(n2)' or var == 'Z2^n2':
                continue

            elif var.startswith('Z3^(n3+'):
                z_power = var.split('+')[1][:-1]
                if z_power == '1':
                    new_z3_term = 'Z3'
                else:
                    new_z3_term = f'Z3^{z_power}'
                new_vars.append(new_z3_term)
            elif var.startswith('Z3^(n3-'):
                z_power = var.split('-')[1][:-1]
                new_z3_term = f'Z3^(-{z_power})'
                new_vars.append(new_z3_term)
            elif var == 'Z3^(n3)' or var == 'Z3^n3':
                continue
            else:
                new_vars.append(var)
        new_terms.append('*'.join(new_vars))
    return ' + '.join(new_terms)

def replace_ni_with_Zi_D_zi(equation):
    terms = top_level_split(equation)
    new_terms = []
    for term in terms:
        variables = top_level_split(term, delimiter='*')
        n1s = 0
        n2s = 0
        n3s = 0
        new_vars = []
        for var in variables:
            if var == 'n1':
                n1s += 1
            elif var == 'n2':
                n2s += 1
            elif var == 'n3':
                n3s += 1
            else:
                new_vars.append(var)
        for i in range(n1s):
            new_vars = new_vars + ['Z1*D_z1']
        for i in range(n2s):
            new_vars = new_vars + ['Z2*D_z2']
        for i in range(n3s):
            new_vars = new_vars + ['Z3*D_z3']
        new_terms.append('*'.join(new_vars))
    return ' + '.join(new_terms)

simplified_DDIs = [pull_out_Z_factor(replace_ni_with_Zi_D_zi(eq)) for eq in filtered_list]

# see how many derivatives get in each equation - affects simplification method
def count_D_number(eq):
    terms = top_level_split(eq)
    all_n_Ds = []
    for term in terms:
        n_Ds = 0
        for D in term:
            if D == 'D':
                n_Ds += 1
        all_n_Ds.append(n_Ds)
    return max(all_n_Ds)

max_n_derivs = max([count_D_number(simplified_DDIs[i]) for i in range(len(simplified_DDIs))])
# max number of derivatives is 2, so can make simple substitution 

def replace_Z1_Dz1_pattern(expr):
    # This regex pattern looks for the literal string "Z1*D_z1*Z1*D_z1"
    pattern = r'Z1\*D_z1\*Z1\*D_z1'
    # Replace with the desired string.
    replacement = '(Z1*D_z1 + Z1^2*D_z1^2)'
    return re.sub(pattern, replacement, expr)

def replace_Z2_Dz2_pattern(expr):
    # This regex pattern looks for the literal string "Z2*D_z2*Z2*D_z2"
    pattern = r'Z2\*D_z2\*Z2\*D_z2'
    # Replace with the desired string.
    replacement = '(Z2*D_z2 + Z2^2*D_z2^2)'
    return re.sub(pattern, replacement, expr)

def replace_Z3_Dz3_pattern(expr):
    # This regex pattern looks for the literal string "Z3*D_z3*Z3*D_z3"
    pattern = r'Z3\*D_z3\*Z3\*D_z3'
    # Replace with the desired string.
    replacement = '(Z3*D_z3 + Z3^2*D_z3^2)'
    return re.sub(pattern, replacement, expr)

def sub_Z_Dz_pattern(expr):
    return replace_Z1_Dz1_pattern(replace_Z2_Dz2_pattern(replace_Z3_Dz3_pattern(expr)))

def combine_Z_powers(equation):
    terms = top_level_split(equation)
    new_terms = []
    for term in terms:
        variables = top_level_split(term, delimiter='*')
        n_Z1s = 0
        n_Z2s = 0
        n_Z3s = 0
        new_vars = []
        D_vars = []
        for var in variables:
            if var == 'Z1':
                n_Z1s += 1
            elif 'Z1^' in var:
                base, power = var.split('^')
                power = Fraction(power.strip('()'))
                n_Z1s += power
            elif var == 'Z2':
                n_Z2s += 1
            elif 'Z2^' in var:
                base, power = var.split('^')
                power = Fraction(power.strip('()'))
                n_Z2s += power
            elif var == 'Z3':
                n_Z3s += 1
            elif 'Z3^' in var:
                base, power = var.split('^')
                power = Fraction(power.strip('()'))
                n_Z3s += power
            elif var.startswith('D'):
                D_vars.append(var)
            else:
                new_vars.append(var)
        if n_Z1s == 1:
            new_vars.append('Z1')
        elif n_Z1s != 0:
            new_vars.append(f'Z1^{n_Z1s}')
        if n_Z2s == 1:
            new_vars.append('Z2')
        elif n_Z2s != 0:
            new_vars.append(f'Z2^{n_Z2s}')
        if n_Z3s == 1:
            new_vars.append('Z3')
        elif n_Z3s != 0:
            new_vars.append(f'Z3^{n_Z3s}')
        new_terms.append('*'.join(new_vars + D_vars))
    return ' + '.join(new_terms)

def write_Z_powers_out(equation):
    terms = top_level_split(equation)
    resulting_terms = []
    for term in terms:
        sub_terms_list = term.split('*')
        for var in sub_terms_list:
            if (var.startswith('Z') or var.startswith('-Z')) and '^' in var:
                Z_var, Z_power = var.split('^')
                if Z_var.startswith('-'):
                    Z_var = Z_var[1:]
                    minus = True
                else:
                    minus = False
                new_Z_string = []
                for i in range(int(Z_power)):
                    new_Z_string.append(Z_var)
                if minus:
                    new_Z_string = '-' + '*'.join(new_Z_string)
                else:
                    new_Z_string = '*'.join(new_Z_string)
                sub_terms_list[sub_terms_list.index(var)] = new_Z_string
        resulting_terms.append('*'.join(sub_terms_list))
    return ' + '.join(resulting_terms)

def write_all_powers_out(equation):
    return pull_minus_signs_to_front(write_Y_powers_out(write_Z_powers_out(equation)))

simplified_DDIs = [fully_expand_equation(sub_Z_Dz_pattern(eq)) for eq in simplified_DDIs]
