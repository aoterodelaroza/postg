#! /usr/bin/octave -q

mol = mol_readxyz("co2_n2o.xyz");
rep = mol_ball(mol);
rep = mol_stick(mol,rep);
rep = mol_stick(mol,rep,"C","N",[3.0 4.0],:,0.02,[0 255 0]);
rep = rep_setdefaultscene(rep,op_rotz(-48) * [eye(3); 0 0 -10]);
rep_write_pov(rep,"co2_n2o.pov");
system("povray -D -UV +Ico2_n2o.pov +Oco2_n2o.png +W1000 +H1000 +A");
