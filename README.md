# Viral Polyprotein Annotation

The core idea for the prediction is the propagation of cleavage sites from already annotated polyproteins to unannotated ones. Annotated polyproteins are first identified in the viral RefSeq database. All viral proteins are then extracted from the RefSeq database and clustered into homologous groups. The groups containing annotated and unannotated polyproteins are aligned. The quality of the alignment and the consistency of cleavage site annotations are evaluated automatically. Finally, the cleavage sites are propagated from annotated to unannotated proteins, and the confidence value of the predictions are determined.


![Alt text](etc/method_workflow.png?raw=true "Title")
