# Wrapper classes to PAML
from __future__ import absolute_import

import os
import re

from phylogenetics.utils import run_subprocess
from dendropy import Tree

class ControlFile:

    def __init__(self, **kwargs):
        """ Base class for writing a PAML control file"""
        # Set the key word arguments as attributes of subclass
        for kw in kwargs:
            setattr(self, kw, kwargs[kw])

        self.input = kwargs

    @property
    def string(self):
        """
            Write control file as single string.
        """
        output_string = ""
        for kw in self.input:
            output_string += kw + " = " + str(self.input[kw]) + "\n"

        return output_string

    def read(self, fname):
        """ Read from outside control file. """
        pass

    def write(self, fname):
        """ Write PAML control file. """
        with open(fname, "w") as f:
            f.write(self.string)


class BaseML(ControlFile):
    """
    Construct a BaseML Control File object
    """
    def __init__(self,
        seqfile="baseml.nuc",
        outfile="baseml.mlb",
        treefile="tree.nwk",
        noisy=0,
        verbose=0,
        runmode=0,
        model=0,
        Mgene=0,
        ndata=1,
        clock=0,
        TipDate="0 100",
        fix_kappa=0,
        kappa=0,
        fix_alpha=0,
        alpha=0,
        Malpha=0,
        ncatG=5,
        fix_rho=1,
        rho=0,
        nparK=0,
        nhomo=0,
        getSE=0,
        RateAncestor=0,
        Small_Diff=1e-6,
        cleandata=1,
        icode=0,
        fix_blength=0,
        method=0
        ):
        """"""
        # Inherit base class
        super(BaseML, self).__init__(
            seqfile=seqfile,
            outfile=outfile,
            treefile=treefile,
            noisy=noisy,
            verbose=verbose,
            runmode=runmode,
            model=model,
            Mgene=Mgene,
            ndata=ndata,
            clock=clock,
            TipDate=TipDate,
            fix_kappa=fix_kappa,
            kappa=kappa,
            fix_alpha=fix_alpha,
            alpha=alpha,
            Malpha=Malpha,
            ncatG=ncatG,
            fix_rho=fix_rho,
            rho=rho,
            nparK=nparK,
            nhomo=nhomo,
            getSE=getSE,
            RateAncestor=RateAncestor,
            Small_Diff=Small_Diff,
            cleandata=cleandata,
            icode=icode,
            fix_blength=fix_blength,
            method=method
        )


    def run(self, ):
        """ Run BaseML """
        pass




class CodeML(ControlFile):
    """
    Construct a CodeML Control File object
    """
    def __init__(self,
        seqfile="sequences.fasta",
        outfile="codeml_output.txt",
        treefile="tree.nwk",
        noisy=3,
        verbose=9,
        runmode=0,
        seqtype=2,
        #CodonFreq=0,
        #ndata=1,
        #clock=0,
        #aaDist=0,
        aaRatefile="lg.dat",
        model=3,
        #NSsites=0,
        #icode=0,
        #Mgene=0,
        #fix_kappa=0,
        #kappa=0,
        #fix_omega=0,
        #omega=0.4,
        fix_alpha=0,
        alpha=1,
        #Malpha=0,
        ncatG=4,
        #fix_rho=1,
        #rho=0,
        #getSE=0,
        RateAncestor=2,
        Small_Diff=1e-6,
        cleandata=0,
        fix_blength=1,
        method=1
        ):
        """
        """
        # Construct path to data files from location of installation
        self.aaRatefile_path = os.path.join(os.path.split(__file__)[0], "dat")

        # Inherit base class
        super(CodeML, self).__init__(
            seqfile=seqfile,
            outfile=outfile,
            treefile=treefile,
            noisy=noisy,
            verbose=verbose,
            runmode=runmode,
            seqtype=seqtype,
            #CodonFreq=CodonFreq,
            #ndata=ndata,
            #clock=clock,
            #aaDist=aaDist,
            aaRatefile= os.path.join(self.aaRatefile_path, aaRatefile),
            model=model,
            #NSsites=NSsites,
            #icode=icode,
            #Mgene=Mgene,
            #fix_kappa=fix_kappa,
            #kappa=kappa,
            #fix_omega=fix_omega,
            #omega=omega,
            fix_alpha=fix_alpha,
            alpha=alpha,
            #Malpha=Malpha,
            ncatG=ncatG,
            #fix_rho=fix_rho,
            #rho=rho,
            #getSE=getSE,
            RateAncestor=RateAncestor,
            Small_Diff=Small_Diff,
            cleandata=cleandata,
            fix_blength=fix_blength,
            method=method
        )



    def run(self, fname="codeml-input.txt"):
        """ """
        # CHECK FOR NECESSARY FILES (PHYLIP and NEWICK)
        if os.path.isfile(self.seqfile) is False:
            raise Exception(""" seqfile does not exist! """)
        if os.path.isfile(self.treefile) is False:
            raise Exception(""" treefile does not exist! """)

        # Write out the control file
        self.write(fname)

        # Run a subprocess for codeml
        run_subprocess("codeml", fname)
