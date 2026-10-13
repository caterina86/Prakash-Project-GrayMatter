function f = get_roi(sourcesubject)

fsdirFROM=sprintf('%s/%s/atlasmgz',cvnpath('freesurfer'),sourcesubject);

% specify the atlas you want to use
surf1file = sprintf('%s/%s.glasser16_atlas.v1_0.mgz%s',fsdirFROM,'lh');
roiL = gifti(surf1file);

surf2file = sprintf('%s/%s.glasser16_atlas.v1_0.mgz%s',fsdirFROM,'rh');
roiR = gifti(surf2file);

f = [squeeze(roiL.cdata)'; squeeze(roiR.cdata)']'; % return as a row vector

end