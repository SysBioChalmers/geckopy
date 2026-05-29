# OpenKineticsPredictor (OKP) direct-API integration

MATLAB GECKO first (active). geckopy port deferred — see the end.

## Goal

Replace the manual OKP workflow (write `OKP.csv` -> upload by hand at
predictor.openkinetics.org -> download `job-*-output.csv` -> parse)
with two functions that talk to the OKP REST API directly:

- `submitOpenKineticsPredictor` — build the input CSV and submit the
  job; persist the returned job id.
- `fetchOpenKineticsPredictor` — check status; if finished, download +
  save + parse into a kcatList. Has a stored-result mode and a polling
  mode.

The existing `writeOpenKineticsPredictorInput.m` and
`readOpenKineticsPredictorOutput.m` are **replaced**: their bodies move
in as local subfunctions of submit (build input) and fetch (parse
output) respectively. The stray `writeOpenKineticsPredictorInput.asv`
is deleted.

## Verified API reference (live, 2026-05-20)

Base: `https://predictor.openkinetics.org/api/v1`
Auth: `Authorization: Bearer ak_...`. Missing/invalid -> HTTP 401
`{"error":"Authentication required. Include the header: ..."}`.
Async model: submit -> poll status -> download result.

### API keys (anonymous, IP-bound)

- `GET /api/api-key/` -> `{"hasKey":bool, "keySuffix":"e86abd"}`
- `POST /api/api-key/generate/` -> `{"key":"ak_...", "keySuffix":"..."}`
  or `{"error":"An active API key already exists for your IP. Revoke it first."}`
- `POST /api/api-key/revoke/` -> `{"revoked":true}`

No login/registration. Key shown only once. Out of scope for the two
functions — the user generates a key in the web UI and provides it.

### POST /submit/ (multipart form)

Form fields: `file=@OKP.csv`, `targets='["kcat"]'`,
`methods='{"kcat":"CataPro"}'`, `handleLongSequences=truncate`,
`includeSimilarityColumns=true`, `canonicalizeSubstrates=true`.

Response:
```json
{"jobId":"BDUASPl","status":"Pending",
 "statusUrl":"/api/v1/status/BDUASPl/",
 "resultUrl":"/api/v1/result/BDUASPl/",
 "quota":{"limit":20000,"used":1,"remaining":19999,"resetsInSeconds":18447}}
```

### GET /status/{jobId}/

```json
{"jobId":"BDUASPl","status":"Completed",
 "submittedAt":"2026-05-20T18:52:32Z","elapsedSeconds":50,
 "queueSeconds":0,"computeSeconds":50,"queuePosition":null,
 "progress":{"moleculesTotal":1,"moleculesProcessed":1,
   "predictionsTotal":1,"predictionsMade":1,"invalidRows":0,
   "stages":[...]},
 "completedAt":"...","resultUrl":"/api/v1/result/BDUASPl/"}
```
- `status` is **Title-case**: `Pending` / `Completed` (also expect
  `Running` / `Failed`). Compare case-insensitively.
- Bad id -> HTTP 404 `{"error":"No job found with id 'NOPE123'."}`.

### GET /result/{jobId}/

- Default: `Content-Type: text/csv`,
  `Content-Disposition: attachment; filename="webkinpred-<jobid>-results.csv"`.
  7 columns:
  `kcat (1/s), Source kcat, Extra Info kcat, Protein Sequence, Substrate, mean similarity to CataPro training data, max similarity to CataPro training data`
- `?format=json`: `{jobId, columns:[...], rowCount, data:[{...}]}` —
  **but** the JSON contains bare `NaN` (`"Extra Info kcat": NaN`),
  which MATLAB `jsondecode` rejects. **Use the CSV form.**

The CSV's columns 1,2,4,5 (kcat, source, sequence, substrate) are
exactly what the existing parser reads; the 2 trailing similarity
columns are ignored by its fixed-index logic.

### POST /validate/ (multipart, optional)

`file=@OKP.csv`, `runSimilarity=false|true` (`true` runs a slow
MMseqs2 pass). Response:
```json
{"rowCount":1,"invalidSubstrates":[],"invalidProteins":[],
 "lengthViolations":{"CataPro":0,...},
 "lengthLimits":{"CataPro":1000,"CatPred":2048,"DLKcat":null,...},
 "similarity":null}
```

### GET /methods/ (no auth) and /quota/ (auth)

`/methods/` lists predictors per target with `id`, `supports`,
`maxSeqLen`, `requiredColumns`. `/quota/` ->
`{limit:20000, used, remaining, resetsInSeconds}` (resets midnight UTC;
failed/partial jobs are credited back).

Live method roster + length limits: CataPro(1000), CatPred(2048),
DLKcat(no limit), EITLEM(1024), KinForm-H/L(1500), MMISA-KM(500),
OmniESI(1000), TurNup(1024), UniKP(1000).

## MATLAB HTTP approach

- **submit** (multipart file upload): `matlab.net.http.RequestMessage`
  with `matlab.net.http.io.MultipartFormProvider` (file part via
  `FileProvider`, the rest as string parts). `webwrite` cannot do
  multipart form uploads, so this is required. Auth via
  `matlab.net.http.HeaderField('Authorization',['Bearer ' apiKey])`.
- **status**: `webread` with `weboptions('HeaderFields',{'Authorization' ['Bearer ' apiKey]})`,
  parse the returned struct.
- **result download**: `websave` (or `webread` returning text) with the
  same auth header; write bytes to `data/OKP_output.csv`.

## Function specs

### submitOpenKineticsPredictor

```matlab
jobId = submitOpenKineticsPredictor(model, ecRxns, modelAdapter, ...
                                    method, apiKey, overwrite)
```
1. Resolve params: `method` arg -> `params.okp.method` -> `'CataPro'`;
   likewise `targets`, `handleLongSequences`, `includeSimilarityColumns`,
   `canonicalizeSubstrates` from `params.okp.*` with built-in defaults.
2. Resolve API key: `apiKey` arg -> env `OKP_API_KEY` ->
   `data/okpApiKey.txt`. Error (naming all three) if none.
3. Build `data/OKP.csv` via the in-lined input builder (former
   `writeOpenKineticsPredictorInput`), honouring `overwrite`.
4. POST `/submit/` (multipart). On non-2xx, error with the `error`
   field.
5. Write `data/OKP_job.txt` (plain text, one field per line):
   ```
   jobId: BDUASPl
   method: CataPro
   targets: kcat
   submittedAt: 2026-05-20T18:52:32Z
   ```
6. Return `jobId`; print the status URL and a "run
   fetchOpenKineticsPredictor later" hint.

### fetchOpenKineticsPredictor

```matlab
[done, kcatList] = fetchOpenKineticsPredictor(model, useStored, ...
                                 modelAdapter, jobId, wait, pollInterval)
```
`useStored` is the **2nd** argument by request.

- `useStored=true`: skip the API entirely. Parse `data/OKP_output.csv`
  (former `readOpenKineticsPredictorOutput` logic) -> `kcatList`,
  `done=true`. Error if the file is absent.
- `useStored=false` (default):
  1. `jobId` arg -> read from `data/OKP_job.txt`.
  2. Resolve API key (same as submit).
  3. GET `/status/{jobId}/`.
     - `Completed` -> GET `/result/{jobId}/` (CSV), save to
       `data/OKP_output.csv`, parse -> `kcatList`, `done=true`.
     - `Pending`/`Running` -> if `wait=true`, poll every `pollInterval`
       seconds until terminal/`timeout`; else print
       `"Job not finished (status: Running). Try later, or check
       https://predictor.openkinetics.org/"`, return `done=false`,
       `kcatList=[]`.
     - `Failed` -> error with the status body.

Returned `kcatList` has the same shape the old reader produced
(`source`, `rxns`, `genes`, `substrates`, `kcats`, `kcatSource` with
`OKP-<method>` labels), consumable by `selectKcatValue`.

## ModelAdapter additions (adapterTemplate.m)

A new optional `okp` block (non-secret config only; the key never goes
here):
```matlab
obj.params.okp.method                   = 'CataPro';
obj.params.okp.targets                  = {'kcat'};
obj.params.okp.handleLongSequences      = 'truncate';
obj.params.okp.includeSimilarityColumns = true;
obj.params.okp.canonicalizeSubstrates   = true;
```
Functions treat `params.okp` as optional and fall back to built-in
defaults, so existing adapters keep working untouched.

## Files

- New: `submitOpenKineticsPredictor.m`, `fetchOpenKineticsPredictor.m`
- Deleted: `writeOpenKineticsPredictorInput.m`,
  `readOpenKineticsPredictorOutput.m`, `writeOpenKineticsPredictorInput.asv`
- Modified: `adapterTemplate.m` (+okp block)

## Decisions locked

- Default method: **CataPro**.
- Load-vs-fetch: **Option A** — one `fetch`, explicit `useStored` (2nd
  arg), no silent fallback. Parsing always happens; fetch returns a
  kcatList.
- API key: arg -> env `OKP_API_KEY` -> `data/okpApiKey.txt`, never the
  adapter.
- Metadata file: plain text `data/OKP_job.txt`, one field per line.
- Result saved to `data/OKP_output.csv`.

## geckopy port (later)

Mirror as `gather_kcats/open_kinetics_predictor/` with `client.py`
(injectable requests.Session), `submit`/`fetch` entry points, the same
build-input / parse-output logic, key via env `OPENKINETICS_API_KEY` /
`.env`. Tests mock HTTP; opt-in live marker. Only needs `requests`
(already a dep). Not started.
