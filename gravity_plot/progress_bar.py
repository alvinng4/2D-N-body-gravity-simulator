import time

import rich.progress


class Progress_bar(rich.progress.Progress):
    def __init__(self):
        super().__init__(
            rich.progress.BarColumn(),
            rich.progress.TextColumn("[green]{task.percentage:>3.0f}%"),
            rich.progress.TextColumn("•"),
            rich.progress.TimeElapsedColumn(),
            rich.progress.TextColumn("•"),
            rich.progress.TimeRemainingColumn(),
        )


class Progress_bar_with_data_size(rich.progress.Progress):
    def __init__(self):
        super().__init__(
            rich.progress.BarColumn(),
            rich.progress.TextColumn("[green]{task.percentage:>3.0f}%"),
            rich.progress.TextColumn("•"),
            rich.progress.TimeElapsedColumn(),
            rich.progress.TextColumn("•"),
            rich.progress.TimeRemainingColumn(),
            "• [magenta]Data size: {task.fields[store_count]}",
        )


def progress_bar_c_lib_adaptive_integrator(tf, t, store_count):
    progress_bar = Progress_bar_with_data_size()
    with progress_bar:
        task = progress_bar.add_task("", total=tf, store_count=1)

        while True:
            # Update progress bar
            progress_bar.update(
                task, completed=t.value, store_count=store_count.value + 1
            )
            time.sleep(0.1)

            if t.value >= tf:
                break

        progress_bar.update(task, completed=tf, store_count=store_count.value + 1)


def progress_bar_c_lib_fixed_integrator(store_npts, store_count):
    progress_bar = Progress_bar_with_data_size()
    with progress_bar:
        task = progress_bar.add_task("", total=store_npts, store_count=1)

        while True:
            # Update progress bar
            progress_bar.update(
                task,
                completed=(store_count.value + 1),
                store_count=(store_count.value + 1),
            )
            time.sleep(0.1)

            if (store_count.value + 1) >= store_npts:
                break

        progress_bar.update(
            task, completed=store_npts, store_count=store_count.value + 1
        )
