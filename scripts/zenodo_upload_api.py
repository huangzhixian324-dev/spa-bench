"""Upload the SPATBench archive to Zenodo via the REST API and publish it.
Usage:  python zenodo_upload_api.py <ACCESS_TOKEN>
Token: Zenodo -> avatar -> Settings -> Applications -> Personal access tokens
       -> New token, scopes: deposit:write, deposit:actions
"""
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ZIP = REPO / 'SPATBench_v33_v34_zenodo.zip'
API = 'https://zenodo.org/api'

TITLE = ('SPATBench: Methodological Choices Can Outweigh Algorithmic Choices in '
         'Transcriptomic Immunotherapy Prediction Benchmarking at Current Sample Sizes')
DESCRIPTION = """<p>SPATBench is an open-source benchmark framework and analysis archive for
transcriptomic immunotherapy-prediction benchmarking.</p>
<p><b>Contents.</b> Complete analysis pipeline and machine-readable results
(36-cell permutation matrix, the 5,000-shuffle Gide 2019 extension, effect-size
bootstrap, one-sample power analysis, XGBoost exclusion evidence), preprocessed
data matrices and cohort metadata, all figures at 300&nbsp;DPI, the manuscript and
supplementary material, and the Circularity Detection Score (CDS) package
v1.2.1.</p>
<p><b>Key results.</b> Switching from a molecular surrogate (cytolytic activity) to
a clinical endpoint (RECIST v1.1) on the same patients collapses every method whose
features overlap the endpoint-defining genes (ElasticNet-Var 0.985 to 0.465; GEP
0.968 to 0.503; IMPRES unchanged), an effect 3.5x the largest method gap.
No trainable method survives FDR correction on clinical endpoints at n&lt;=43,
whereas at n=73 (Gide 2019) ElasticNet-MI is significant (p=0.0232, 5,000 shuffles).</p>
<p><b>Code.</b> https://github.com/huangzhixian324-dev/spa-bench and
https://github.com/huangzhixian324-dev/cds (MIT).</p>"""

METADATA = {
    'metadata': {
        'title': TITLE,
        'upload_type': 'software',
        'publication_date': time.strftime('%Y-%m-%d'),
        'creators': [{'name': 'Huang, Zhixian',
                      'affiliation': 'Putian University'}],
        'description': DESCRIPTION,
        'access_right': 'open',
        'license': 'MIT',
        'keywords': ['immunotherapy', 'benchmarking', 'circularity',
                     'endpoint definition', 'transcriptomics',
                     'permutation testing', 'reproducibility'],
        'related_identifiers': [
            {'relation': 'isSupplementTo',
             'identifier': 'https://github.com/huangzhixian324-dev/spa-bench',
             'resource_type': 'software'},
        ],
    }
}


def req(method, url, token, data=None, headers=None, raw=None):
    h = {'Authorization': f'Bearer {token}'}
    if headers:
        h.update(headers)
    body = None
    if raw is not None:
        body = raw
    elif data is not None:
        body = json.dumps(data).encode()
        h['Content-Type'] = 'application/json'
    r = urllib.request.Request(url, data=body, headers=h, method=method)
    try:
        with urllib.request.urlopen(r, timeout=600) as resp:
            text = resp.read().decode('utf-8', errors='replace')
            return resp.status, (json.loads(text) if text.strip().startswith(('{', '[')) else text)
    except urllib.error.HTTPError as e:
        return e.code, e.read().decode('utf-8', errors='replace')[:400]


def main(token):
    if not ZIP.exists():
        print('zip 不存在:', ZIP); return
    print(f'上传包: {ZIP.name}  {ZIP.stat().st_size/1024/1024:.1f} MB\n')

    # 1) create deposition
    st, dep = req('POST', f'{API}/deposit/depositions', token, data={})
    print('[1] create deposition:', st)
    if st not in (200, 201):
        print('   ', dep); return
    dep_id = dep['id']
    bucket = dep.get('links', {}).get('bucket')
    print(f'    deposition id = {dep_id}')

    # 2) upload file (bucket API, supports large files)
    if bucket:
        upload_url = f'{bucket}/{ZIP.name}'
    else:
        upload_url = f'{API}/deposit/depositions/{dep_id}/files'
    print('[2] uploading...')
    with open(ZIP, 'rb') as f:
        st, out = req('PUT', upload_url, token, raw=f.read(),
                      headers={'Content-Type': 'application/octet-stream'})
    print('    upload status:', st)
    if st not in (200, 201):
        print('   ', str(out)[:300]); return

    # 3) metadata
    st, out = req('PUT', f'{API}/deposit/depositions/{dep_id}', token, data=METADATA)
    print('[3] metadata:', st)
    if st not in (200, 201):
        print('   ', str(out)[:300])

    # 4) publish
    st, out = req('POST', f'{API}/deposit/depositions/{dep_id}/actions/publish', token)
    print('[4] publish:', st)
    if st in (200, 202):
        doi = out.get('doi') if isinstance(out, dict) else None
        print('\n✅ 已发布！DOI =', doi)
        print('   URL =', out.get('links', {}).get('record_html') if isinstance(out, dict) else '')
        (REPO / 'results' / 'zenodo_doi.txt').write_text(str(doi), encoding='utf-8')
    else:
        print('   ', str(out)[:400])
        print('\n（如提示未发布，可登录 Zenodo 在 My dashboard 手动点 Publish；文件已上传）')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('用法: python zenodo_upload_api.py <ACCESS_TOKEN>')
        print('token: Zenodo -> Settings -> Applications -> Personal access tokens')
    else:
        main(sys.argv[1].strip())
