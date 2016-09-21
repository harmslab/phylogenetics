Metadata format
===============

A common issue today when doing phylogenetics
is managing many different data formats that exist (i.e. fasta, phylip,
newick, nexus, paml, etc.) while the data remains persistent. `phylogenetics` deals with this issue by tranforming
all data that is read/written by the package to a single defined format metadata
format. It follows a simple key-value dictionary structure. This data format is
lightweight and easy to dump as pickle, json, csv, etc. This also allows us to
continue developing the API without fear of backward incompatibility.

Check out the different metadata formats:

* `Project metadata format`_
* `HomologSet metadata`_
* `Alignment metadata format`_
* `Tree metadata format`_
* `AncestorSet metadata format`_


Project metadata format
-----------------------
::

    {
        "name" : "human-project",
        "date" : "",
        "edited" : "",
        "SequenceList" :
            {
                "history" : [],
                "contents" : [
                    {
                        "id" : "XX00000000",
                        "sequence" : "AGAMAMGATKLLSMA",
                        "orgname" : "human",
                    },
                    {
                        "id" : "XX00000001",
                        "sequence" : "AGAKKLGATKLLSMA",
                        "orgname" : "human",
                    }
                ]
            },
        "AlignmentList" : {
            "id" : "AlignmentList0",
            "history" : [],
            "contents" : [
                {
                    "id" : "Alignment0",
                    "history" : [],
                    "contents" : [
                        {
                            "id" : "Ali0000000",
                            "sequence" : "AGAMAMGATKLLSMA",
                            "Sequence" : "Seq0000000"
                        },
                        {
                            "id" : "Ali0000001",
                            "sequence" : AGAKKLGATKLLSMA"",
                            "Sequence" : "Seq0000001"
                        },
                    ],
                },
                {
                    "name" : "old",
                    "history" : [],
                    "alignment" : [
                        {
                            "id" : "Ali0000000",
                            "sequence" : "",
                            "Sequence" : "Seq0000000"

                        },
                        {
                            "id" : "Ali0000001",
                            "sequence" : "",
                            "Sequence" : "Seq0000001"

                        },
                    ],
                }
            ],
        },
        "TreeList" : {
            "id" : "TreeList0"
            "contents" : [
                {
                    "id" : "tree000000",
                    "date" : "",
                    "notes" : "This is the best tree.",
                    "Alignment" : "Alignment0",
                    "stats" : {
                        "supports" : "aLRT",
                    }
                    "newick" : ((,),),
                },
                {
                    "name" : "tree1",
                    "date" : "",
                    "notes" : "This is the old tree.",
                    "stats" : {
                        "supports" : "SH",
                    }
                    "newick" : ((,),),
                }
            ],
        },
        "AncestorList" : [
            {
                "tree" : "tree0",
                "notes" : "",
                "date" : "",
                "ancestors": [
                    {
                        "id" : "anc0000000",
                        "mlsequence" : "AGAMAMGATKLLSMA",
                        "posterior" : [],
                    },
                    {
                        "id" : "anc0000001",
                        "mlsequence" : "AGAKKLGATKLLSMA",
                        "posterior" : [],
                ],
            },
            {
                "tree" : "tree0",
                "notes" : "",
                "date" : "",
                "ancestors": [
                    {
                        "id" : "anc0000000",
                        "mlsequence" : "AGAMAMGATKLLSMA",
                        "posterior" : [],
                    },
                    {
                        "id" : "anc0000001",
                        "mlsequence" : "AGAKKLGATKLLSMA",
                        "posterior" : [],
                ],
            },
        ],
    }

HomologSet metadata
-------------------
::

    {
        "name" : "dataset1",
        "date" : "",
        "edited" : "",
        "homologs" : [
            {
                "id" : "XX00000000",
                "sequence" : "AGAMAMGATKLLSMA",
                "orgname" : "human",
            },
            {
                "id" : "XX00000001",
                "sequence" : "AGAKKLGATKLLSMA",
                "orgname" : "human",
            }
        ]
    }

Alignment metadata format
-------------------------
::

    {
        "name" : "align0",
        "latest": True,
        "date" : "",
        "alignment" : [
            {
                "id" : "seq0000000",
                "sequence" : "",
            },
            {
                "id" : "seq0000001",
                "sequence" : "",
            },
        ],
    }


Tree metadata format
----------------------
::

    {
        "name" : "tree0",
        "date" : "",
        "alignment" : "align0",
        "notes" : "This is the best tree.",
        "stats" : {
            "supports" : "aLRT",
        }
        "newick" : ((,),),
    }


AncestorSet metadata format
---------------------------
::

    {
        "tree" : "tree0",
        "notes" : "",
        "date" : "",
        "ancestors": [
            {
                "id" : "anc0000000",
                "mlsequence" : "AGAMAMGATKLLSMA",
                "posterior" : [],
            },
            {
                "id" : "anc0000001",
                "mlsequence" : "AGAKKLGATKLLSMA",
                "posterior" : [],
        ],
    },
