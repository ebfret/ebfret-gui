function close(self)
    % get rid of ui elements
    close(self.handles.mainWindow);
    % get rid of object data structure
    delete(self);