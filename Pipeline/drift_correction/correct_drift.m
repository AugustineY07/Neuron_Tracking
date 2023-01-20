function correct_drift(z_drift,locdir2,id,ish)
        centXZY = readNPY([locdir2,'\centXZY.npy']);
        centXZY_corrected = centXZY;
        centXZY_corrected(:,2) = centXZY(:,2) - z_drift;
        writeNPY(centXZY_corrected,[locdir2,'\centXZY_corrected.npy'])
end


