/*******************************************************************************
 * Copyright (C) 2022-2023 Simone Rubinacci
 * Copyright (C) 2022-2023 Olivier Delaneau
 *
 * MIT Licence
 ******************************************************************************/

#ifndef _RETRY_IO_H
#define _RETRY_IO_H

#include <string>
#include <chrono>
#include <thread>

#include <utils/otools.h>

//Outcome of one attempt passed to retry_with_backoff.
struct attempt_result
{
	bool ok = false;      //attempt succeeded: stop and return
	bool fatal = false;   //failure is non-retryable: abort immediately, do not retry
	std::string err;      //human-readable error, used in the warning/error messages
};

//Retry `attempt` up to n_retry times with exponential backoff. With n_retry=3 and
//base_delay=1s the attempts are separated by 1s then 2s sleeps (no sleep before the
//first attempt). Between failed attempts a warning is logged with the correct number of
//attempts still remaining. A fatal outcome aborts immediately via vrb.error; exhausting
//all attempts also aborts via vrb.error. `what` names the action for the log messages,
//e.g. "reading binary reference panel [path]".
//
//Note: this is a new phase-local header. If ligate/concordance ever grow the same
//retry need, we should move it to common/src/utils/ with symlinks
    — but I kept it phase-local to avoid touching the shared tree unnecessarily.
template <typename Fn>
void retry_with_backoff(const std::string& what, int n_retry, std::chrono::seconds base_delay, Fn&& attempt)
{
	std::chrono::seconds delay = base_delay;
	for (int i = 0; i < n_retry; ++i)
	{
		const attempt_result r = attempt();
		if (r.ok) return;
		if (r.fatal) vrb.error("Non-retryable error while " + what + ": " + r.err);

		const int attempts_left = n_retry - i - 1;
		if (attempts_left > 0)
		{
			vrb.warning("Error while " + what + ": " + r.err + ". " + std::to_string(attempts_left) + " attempt(s) remaining");
			std::this_thread::sleep_for(delay);
			delay *= 2;
		}
		else vrb.error("Max number of retries attempted while " + what + ". Last error: " + r.err);
	}
}

#endif
