# Codex Task Brief Template (Non‑Programmer Friendly)

Use this template at the start of each new Codex chat.
Replace the placeholders in `[brackets]`.

---

## 1) Task Goal (plain language)
I want you to do exactly this:
- [one concrete outcome]
- [optional second outcome]

Success means:
- [what I can see/hear/verify when done]

## 2) Scope (files/folders allowed)
You may edit only:
- [file/folder A]
- [file/folder B]

Do NOT edit:
- [protected file/folder A]
- [protected file/folder B]

## 3) Constraints (must follow)
- Keep changes minimal and focused.
- Preserve backward compatibility for presets/snapshots/parameter IDs unless I explicitly approve changes.
- Keep audio thread RT-safe (no locks, no allocations, no blocking).
- If this task affects latency/timing/routing, call it out explicitly.

## 4) Deliverables
Provide:
1. Short plan before editing.
2. Code/docs patch.
3. Short “Before vs After” explanation in non-technical language.
4. Risks and rollback note.

## 5) Validation commands
Run these checks and report results:
- [build command]
- [test command]
- [targeted grep/sed command]

If no automated tests exist, run the best available sanity checks and say so clearly.


## 5b) Build policy in Codex/cloud (environment-specific)
Use this policy when build steps may fetch JUCE from GitHub:
- Run preflight first: `cmake -S . -B build`
- If output shows JUCE FetchContent clone/fetch HTTP 403 (or CONNECT tunnel failed):
  - mark build as `⚠️ blocked by known network restriction (JUCE FetchContent GitHub access denied in Codex/cloud)`
  - do not retry full build in the same run
  - continue with local-only validations (grep/sed/unit tests that require no network)
  - in final report include the exact failing command + first relevant error lines

## 6) Output format I want
At the end, respond with:
- Summary (bullet points)
- Files changed
- Testing results (✅/⚠️/❌ per command)
- Open questions for me (if any)

## 7) Safety gate (important)
If you discover this task is larger than expected, STOP and give me:
- a reduced safe patch option,
- a full patch option,
- and your recommendation.

---

## Optional add-ons you can append per task
- Audio behavior focus: [e.g., preserve transients / avoid pumping]
- Standards focus: [EBU R128 / A/85]
- Channel/layout focus: [stereo / 5.1 / 7.1.4]
- Performance budget: [CPU % target]
- Diff budget: [max files or rough line count]


## 8) Conflict Playbook (GitHub merge safety)
Use this playbook whenever a PR/merge conflict appears.

### A) Before requesting a patch
- Work from a branch that is updated from the latest `main`.
- Ask for a **minimal, block-scoped patch** (no full-file rewrites unless explicitly requested).
- In the prompt, specify exact headings/keys that may be edited.

### B) Conflict-resistant patch request pattern
Include these hard rules in prompts:
- "Do NOT rewrite the whole file."
- "Edit only the exact block/keys listed below."
- "Preserve all lines not in scope."
- "Return unified diff for changed file(s) only."

### C) If conflict markers appear in GitHub
Conflict markers look like:
- `<<<<<<<` current branch
- `=======`
- `>>>>>>>` incoming branch

Resolution steps:
1. Keep the desired final content.
2. Delete all three marker lines.
3. Re-check indentation/syntax for the file type (YAML uses `#` for comments, not `//`).
4. Save and run a quick sanity check.

### D) Docs conflict fallback (safe mode)
If repeated conflicts occur in the same doc file:
- Request from Codex the **full final content for that one file only**.
- Replace that file fully once in GitHub editor.
- Avoid mixing unrelated edits in the same merge.

### E) Branch hygiene
- Prefer one short-lived branch per task.
- Merge quickly after review, then delete merged branches.
- Keep one optional `integration/docs` branch only if combining multiple doc PRs.

### F) Mandatory post-merge checks
- Confirm no conflict markers remain: `rg -n "^(<<<<<<<|=======|>>>>>>>)" Docs`
- Validate syntax for edited formats (YAML/JSON/Markdown).
- Record exact validation commands in PR notes.

