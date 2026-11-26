"""Detect accidental filename/content substitutions and missing example dependencies."""
from pathlib import Path
import json,hashlib,sys
ROOT=Path(__file__).resolve().parents[1]
def sha(data):return hashlib.sha256(data).hexdigest()
def statements(text):return '\n'.join(line.rstrip() for line in text.splitlines() if line.strip() and not line.lstrip().startswith('!'))
def main():
    manifest=json.loads((ROOT/'source_manifest.json').read_text(encoding='utf8'))
    expected=set();fail=[]
    for row in manifest['files']:
        p=ROOT/row['path'];expected.add(row['path'])
        if not p.is_file():fail.append('Missing '+row['path']);continue
        raw=p.read_bytes();text=raw.decode(row['original_encoding'])
        if sha(raw)!=row['sha256']:fail.append('Content checksum differs: '+row['path'])
        if sha(raw)!=row['original_sha256'] or p.name!=row['original_name']:fail.append('Original filename/bytes not preserved: '+row['path'])
        if sha(statements(text).encode())!=row['statement_sha256']:fail.append('Numerical statements differ: '+row['path'])
        for needle in row['required_content']:
            if needle.lower() not in text.lower():fail.append('Role/content mismatch: '+row['path'])
    actual={p.relative_to(ROOT).as_posix() for folder in ['src','tests'] for p in (ROOT/folder).rglob('*.f90')}
    if actual!=expected:fail.append('Unlisted or missing Fortran sources: '+str(sorted(actual^expected)))
    if list(ROOT.glob('*.f90')):fail.append('Unexpected top-level legacy Fortran file')
    for row in manifest['fixtures']:
        p=ROOT/row['path']
        if not p.is_file() or sha(p.read_bytes())!=row['sha256']:fail.append('Fixture differs: '+row['path'])
    case=ROOT/'examples/smoke_2_2_2'
    lines=(case/'input.in').read_text().splitlines()
    for name in lines[:6]+[lines[10]]:
        if not (case/name).is_file():fail.append('Missing input dependency: '+name)
    if fail:print('\n'.join(fail));raise SystemExit(1)
    print(f'SOURCE_INTEGRITY_PASSED sources={len(expected)} fixtures={len(manifest["fixtures"])}')
if __name__=='__main__':main()
