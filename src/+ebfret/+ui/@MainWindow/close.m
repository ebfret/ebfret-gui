function close(self)
    % close ui elements
    % (we are assuming matlab deletes dependent ui elements correctly)
    delete(self.handles.mainWindow);
    % delete containing class structure
    delete(self);
