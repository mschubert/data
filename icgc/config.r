version = "release_23"
raw_data = module_file('./data', version, mustWork=TRUE)
cached_data = file.path(module_file('./cache'), version)

if (!file.exists(cached_data))
    dir.create(cached_data)
