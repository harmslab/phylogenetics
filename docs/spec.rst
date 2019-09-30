

.. ::
    project
        metadata.json
        0.csv
        1.csv
        2.csv
        ...


.. ::python

    p = Project(project_dir)
    p.load_fasta("some.fasta")
    p.blast()
    p.align()
    p.build_tree()
    
    project.set_current(3)
    
   
At every step, it says "writing project_dir/junk.csv"
metadata.json
    
if df changed due to manual operations, save out on fly

wrap 
