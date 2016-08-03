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
        "HomologSet" : [
            {
                "id" : "seq00000000",
                "sequence" : "AGAMAMGATKLLSMA",
                "orgname" : "human",
            },
            {
                "id" : "seq00000001",
                "sequence" : "AGAKKLGATKLLSMA",
                "orgname" : "human",
            }
        ],
        "Alignment" : [
            {
                "seq00000000" : "AG---KKLGATKL"
            },
            {
                "seq00000001" : ""
            },
        ],
        ""
    }

HomologSet metadata
-------------------
::

    {
        "HomologSet" : [
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
        "Alignment" : {
            "latest" : [
                {
                    "id" : "seq0000000",
                    "sequence" : "",
                    "date" : "",
                },
                {
                    "id" : "seq0000001",
                    "sequence" : "",
                    "data" : "",
                },
            ],
            "align0" : [
                {
                    "id" : "seq0000000",
                    "sequence" : "",
                },
                {
                    "id" : "seq0000001",
                    "sequence" : "",
                },
            ]
        }
    }


Tree metadata format
--------------------
::

    {
        "Trees" : [
            {
                "name" : "tree0",
                "notes" : "This is the best tree.",
                "stats" : {
                    "supports" : "aLRT",
                }
            },
            {
                "name" : "tree1",
                "notes" : "This is a bad tree.",
                "stats" : {
                    "supports" : "SH"
                }
            }
        ]
    }


AncestorSet metadata format
---------------------------
::

    {
        "Ancestors" : [
            {
                "treename" : "tree0",
                "notes" : "",
                "date" : "",
                "AncestorSet": [
                    {
                        "id" : "anc0000000",
                        "mlsequence" : "AGAMAMGATKLLSMA",
                    },
                    {
                        "id" : "anc0000001",
                        "mlsequence" : "AGAKKLGATKLLSMA",
                ],
            },
        ]
    }
