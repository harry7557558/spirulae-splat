"""Single toggle controlling the per-iteration training profiler.

When True:
  - trainer._train_step inserts torch.cuda.synchronize() probes and prints
    [PROF n=...] every 100 iterations after a 10-iter warmup.
  - model.engine_train_step inserts a pre-call synchronize() and prints
    [PROF_ETS n=...] with the camera/gt/weights/presync/ccall/post split.
  - datamanager.next_train prints [PROF_NT n=...] with the
    pre/get_batch/warp/pack/post split.

When False: zero overhead — no perf_counter_ns calls, no syncs, no prints.

Flip this constant to True for one run, restart training, copy the numbers,
flip back to False.
"""

PROFILE_TRAIN_STEP = False
