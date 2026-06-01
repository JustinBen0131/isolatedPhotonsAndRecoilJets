# Meeting Capture And Zoom

This policy covers Justin's explicit requests to create Zoom rooms, capture
meeting context, and convert meeting output into ThesisAnalysis tasks.

## Trigger Phrase

When Justin says a close variant of:

- "set up a Zoom room";
- "open a Zoom room";
- "make a Zoom room";
- "start a Zoom room";
- "get me a Zoom invite link";

Codex should treat that as explicit permission to use Computer Use on the local
Zoom app or Zoom web UI for the narrow purpose of creating a room and returning
the invite link.

## Zoom Room Setup Contract

For a requested ad hoc Zoom room, Codex should:

1. Open Zoom with Computer Use.
2. Start or create the appropriate Justin-owned meeting room.
3. Enable or verify available transcript and AI-summary capture for that room:
   cloud recording, audio transcript, AI Companion meeting summary, smart
   chapters/highlights, and next steps where Zoom exposes those controls.
4. Copy the invite link or invitation text.
5. Return only the clean invite link and any short caveats Justin needs before
   sending it to the person joining.

If Zoom requires a host/admin setting, payment, account upgrade, credential,
or permission gate, stop and report the exact blocker. Do not ask for or handle
passwords/secrets.

## Privacy And Visibility

Do not imply that transcription is invisible. Zoom recording, transcription,
and AI Companion features can notify participants or be visible in the meeting
UI depending on host/account settings. Assume collaborators may see that a
meeting is recorded or summarized.

For meetings Justin hosts, Codex may enable/verify Justin-owned Zoom capture
settings when Justin asks. For meetings Justin only joins, Codex should not try
to secretly record, transcribe, or summarize. Instead, use available host
artifacts, calendar links, Indico pages, slides, emails, Justin notes, or a
host-approved transcript.

## Post-Meeting Sweep

For meetings with transcript/summary artifacts, Codex should convert the
meeting into a compact OS update:

- extract decisions, action items, open questions, named collaborator comments,
  plot requests, approval risks, and slide/document follow-ups;
- save or verify the meeting slide PDF under
  `/Users/patsfan753/Desktop/ThesisAnalysis/Presentations` with a dated,
  organized name when accessible;
- update `agent_context/CODEX_WORK_REGISTER.yaml` first;
- sync Linear workstream issues second;
- update Today's Plan only for cockpit-worthy active/waiting/review items;
- keep raw transcripts and AI summaries private unless Justin explicitly
  promotes sanitized notes.

If transcript processing is not ready, record a waiting status and next check
instead of treating the artifact as absent.

## External Action Guardrails

Do not email collaborators, invite bots, mutate Google Slides, submit SDCC
work, publish transcripts, or change account-level Zoom/admin settings unless
Justin explicitly asks for that specific action.
