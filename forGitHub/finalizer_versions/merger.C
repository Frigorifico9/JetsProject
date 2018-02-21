{
	TFileMerger merger;
	for (int i = 0 ; i < 500 ; ++i)
	{
		merger->AddFile(Form("/storage/alice/frigorifico/jetsProject/JOBS/%d/histos.root",i));
	}
	merger->OutputFile("mergedHistos.root");
	merger->Merge();

}
