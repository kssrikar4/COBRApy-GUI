import re
import difflib

_PAREN_RE = re.compile(r"\(([^()]+)\)\s*$")

def extract_id_candidate(s: str) -> str:
    """Extract ID candidate from string with parentheses."""
    m = _PAREN_RE.search(s.strip())
    return m.group(1).strip() if m else s.strip()

def resolve_reaction_id(s: str, model, fuzzycutoff: float = 0.6, maxsuggestions: int = 8):
    """
    Always return a valid reaction ID from a user/display string (or raise ValueError with suggestions).
    Resolution order:
      1) Exact ID
      2) Extracted '(id)' from 'name (id)'
      3) Exact name (unique)
      4) Substring over ids+labels (unique)
      5) Fuzzy over ids+labels (unique)
    """
    q = s.strip()
    ids = [r.id for r in model.reactions]
    labels = [f"{r.name} ({r.id})" if r.name else r.id for r in model.reactions]
    label_to_id = {lab: rid for lab, rid in zip(labels, ids)}
    name_to_id = {}
    for r in model.reactions:
        if r.name:
            nm = r.name.strip()
            if nm not in name_to_id:
                name_to_id[nm] = r.id

    if q in model.reactions:
        return q
    qid = extract_id_candidate(q)
    if qid in model.reactions:
        return qid
    if q in name_to_id:
        return name_to_id[q]
    qlow = q.lower()
    subs_ids = [rid for rid in ids if qlow in rid.lower()]
    subs_labs = [lab for lab in labels if qlow in lab.lower()]
    subs = list({*(subs_ids), *[label_to_id[lab] for lab in subs_labs]})
    if len(subs) == 1:
        return subs[0]
    f_ids = difflib.get_close_matches(q, ids, n=max_suggestions, cutoff=fuzzy_cutoff)
    f_labs = difflib.get_close_matches(q, labels, n=max_suggestions, cutoff=fuzzy_cutoff)
    fuzz = list({*f_ids, *[label_to_id[lab] for lab in f_labs]})
    if len(fuzz) == 1:
        return fuzz[0]
    suggestions = (subs or fuzz)[:max_suggestions]
    hint = f" Suggestions: {suggestions}" if suggestions else ""
    raise ValueError(f"Could not resolve a valid reaction ID from: '{s}'.{hint}")

def resolve_ddca_shortcut(s: str, model):
    """Resolve laurate/ddca shortcut to reaction IDs."""
    key = s.strip().lower()
    aliases = {"n-c12:0", "c12:0", "laurate", "dodecanoate", "ddca"}
    if key in aliases:
        picks = [rid for rid in ["EX_ddca_e", "DDCAt2pp"] if rid in model.reactions]
        if picks:
            return picks
    return []

def extract_id(option_string):
    """Extract ID from option string using regex."""
    match = re.search(r'\((.*?)\)', option_string)
    if match:
        return match.group(1)
    return option_string
