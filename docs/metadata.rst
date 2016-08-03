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



Project metadata format
-----------------------
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
        ],
        "Alignment" : [
            {
                "XX00000000" : "AG---KKLGATKL"
            },
            {
                "XX00000001" : ""
            },
        ],
        ""
    }

HomologSet metadata
-------------------
::

    {"HomologSet" : [
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

    {"Alignment" : {
            "latest" : [
                {
                    "id" : "XX00000000"
                },
                {},
            ]

        }
    }
