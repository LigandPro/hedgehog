from hedgehog.stages.docking.utils import run_docking


def main(config: dict, reporter=None):
    return run_docking(config, reporter=reporter)
