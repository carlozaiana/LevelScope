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
