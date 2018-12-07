import click
from click import (
    Path,
)

from ._kmkm import PyKmerCounter as KmerCounter
from .collection import KmerCollection
from .logger import LOGGER as LOG, enable_logging, DEBUG

def handle_logging_args(verbose, quiet):
    if quiet:
        return
    if verbose > 0:
        enable_logging(DEBUG)
    else:
        enable_logging()


@click.group()
def main():
    pass


@click.command("count")
@click.argument('outfile', required=True, type=Path())
@click.argument('seqfiles', nargs=-1, required=True, type=Path(exists=True))
@click.option('-k','--ksize', default=21, type=int)
@click.option('-c', '--cvsize', default=100000000, type=int)
@click.option('-v', '--verbose', count=True)
@click.option('-q', '--quiet', default=False)
def count_file(outfile, seqfiles, ksize, cvsize, quiet, verbose):
    handle_logging_args(verbose, quiet)
    LOG.info("Counting files...")
    kc = KmerCounter(ksize, cvsize)
    for sf in seqfiles:
        LOG.info("\t" + sf)
        kc.count_file(sf)
    LOG.info("Saving to " + outfile)
    kc.save(outfile)
    LOG.info("All done!")
main.add_command(count_file)

@click.command("counteach")
@click.argument('outfile', required=True, type=Path())
@click.argument('seqfiles', nargs=-1, required=True, type=Path(exists=True))
@click.option('-k','--ksize', default=21, type=int)
@click.option('-c', '--cvsize', default=100000000, type=int)
@click.option('-a', '--append', default=False, is_flag=True)
@click.option('-v', '--verbose', count=True)
@click.option('-q', '--quiet', default=False, is_flag=True)
def counteach(outfile, seqfiles, ksize, cvsize, append, quiet, verbose):
    handle_logging_args(verbose, quiet)
    LOG.info("Counting sequence files...")
    mode = "a" if append else "w"
    kc = KmerCollection(outfile, mode=mode, ksize=ksize, cvsize=cvsize)
    for f in seqfiles:
        kc.count_seqfile(f)
    LOG.info("All done!")
main.add_command(counteach)

@click.command("aggregate")
@click.argument('outfile', required=True, type=Path())
@click.argument('countfiles', nargs=-1, required=True, type=Path(exists=True))
@click.option('-a', '--append', default=False, is_flag=True)
@click.option('-v', '--verbose', count=True)
@click.option('-q', '--quiet', default=False, is_flag=True)
def aggregate(outfile, countfiles, append, quiet, verbose):
    handle_logging_args(verbose, quiet)
    LOG.info("Counting files...")
    mode = "a" if append else "w"
    kc = KmerCollection(outfile, mode=mode)
    for cf in countfiles:
        kc.add_file(cf)
    LOG.info("All done!")
main.add_command(aggregate)

