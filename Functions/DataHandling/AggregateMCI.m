MCIAll = struct;

disp("%%%%%%%%%%%%%% Aggregating MCI with the following order MCIPos, MCINeg, MCI Unknown %%%%%%%%%%%%%%");
MCIAll.FlagPos          = [MCIPos.FlagPos MCINeg.FlagPos Unknown.FlagPos];
MCIAll.FlagSpawnTime    = [MCIPos.FlagSpawnTime MCINeg.FlagSpawnTime Unknown.FlagSpawnTime];
MCIAll.FlagTrigPos      = [MCIPos.FlagTrigPos MCINeg.FlagTrigPos Unknown.FlagTrigPos];
MCIAll.FlagTrigTimes    = [MCIPos.FlagTrigTimes MCINeg.FlagTrigTimes Unknown.FlagTrigTimes];
MCIAll.TrigPos          = [MCIPos.TrigPos MCINeg.TrigPos Unknown.TrigPos];
MCIAll.TrigPosTimes     = [MCIPos.TrigPosTimes MCINeg.TrigPosTimes Unknown.TrigPosTimes];
MCIAll.OutOfBoundPos    = [MCIPos.OutOfBoundPos MCINeg.OutOfBoundPos Unknown.OutOfBoundPos];
MCIAll.PathKey          = [MCIPos.PathKey MCINeg.PathKey Unknown.PathKey];
MCIAll.Path             = [MCIPos.Path MCINeg.Path Unknown.Path];
MCIAll.Cond             = [MCIPos.Cond MCINeg.Cond Unknown.Cond];
MCIAll.CondTable        = [MCIPos.CondTable MCINeg.CondTable Unknown.CondTable];
MCIAll.Info             = [MCIPos.Info MCINeg.Info Unknown.Info];
MCIAll.Errors           = [MCIPos.Errors MCINeg.Errors Unknown.Errors];
MCIAll.ReconstructedOOB = [MCIPos.ReconstructedOOB MCINeg.ReconstructedOOB Unknown.ReconstructedOOB];
